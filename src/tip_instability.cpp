/* Fused tip-instability kernel.
 *
 * Replaces, for the batch path of R's TipInstability() (log = TRUE, uniform
 * tree dimensions), the chain
 *     .Call(LOG_GRAPH_GEODESIC_MULTI)  ->  R matrix()  ->  Rfast row-stats
 *     ->  symmetric matrix  ->  Rfast::rowmeans
 * with a single call that returns the length-nTip instability vector directly.
 *
 * The nPairs x nTrees distance matrix is never materialised as an R object:
 *  - mean + sd  : streamed per pair via Welford (O(nPairs) memory);
 *  - median/mad : held only in a transient C++ buffer, reduced with
 *                 std::nth_element (no Rfast dependency).
 *
 * Numerics mirror the previous Rfast path exactly: sd uses the (n-1) divisor,
 * mad uses the 1.4826 constant, medians use the mean-of-two-middle rule for
 * even n, and non-finite deviations collapse to 0.
 */

#include <R.h>
#include <Rinternals.h>
#include <vector>
#include <algorithm>
#include <cmath>
#include "geodesic.h"

/* Median of [first, last) via nth_element; reorders the range in place.
 * Even n: mean of the two central order statistics (matches matrixStats /
 * the Rfast path this replaces). */
static double median_inplace(double *first, double *last) {
  const std::ptrdiff_t n = last - first;
  if (n <= 0) return 0.0;
  const std::ptrdiff_t mid = n / 2;
  std::nth_element(first, first + mid, last);
  const double hi = first[mid];
  if (n & 1) {
    return hi;
  }
  const double lo = *std::max_element(first, first + mid);
  return 0.5 * (lo + hi);
}

extern "C" SEXP TIP_INSTABILITY(SEXP n_tip, SEXP n_node, SEXP parent_all,
                                SEXP child_all, SEXP n_edge, SEXP n_tree,
                                SEXP which_ave, SEXP which_dev) {
  const int nTip = INTEGER(n_tip)[0];
  const int nNode = INTEGER(n_node)[0];
  const int allNodes = nTip + nNode;
  const int nEdge = INTEGER(n_edge)[0];
  const int nTree = INTEGER(n_tree)[0];
  const int whichAve = INTEGER(which_ave)[0];  /* 1 = mean,  2 = median */
  const int whichDev = INTEGER(which_dev)[0];  /* 1 = sd,    2 = mad    */
  const R_xlen_t nPairs = (R_xlen_t)nTip * (nTip - 1) / 2;
  const int *parAll = INTEGER(parent_all);
  const int *chAll = INTEGER(child_all);

  /* Reused across trees by graph_geodesic_phylo(). */
  std::vector<int> interim((size_t)allNodes * allNodes);

  /* Per-pair deviation and average, in column-major lower-triangle order
   * (outer column j, inner row i > j) -- the same ordering R uses. */
  std::vector<double> devs(nPairs), aves(nPairs);

  const bool needBuffer = (whichAve == 2 || whichDev == 2);

  if (!needBuffer) {
    /* mean + sd: stream Welford accumulators, no nPairs x nTrees buffer. */
    std::vector<double> mean(nPairs, 0.0), m2(nPairs, 0.0);
    for (int t = 0; t < nTree; ++t) {
      R_CheckUserInterrupt();
      const int *par = parAll + (R_xlen_t)t * nEdge;
      const int *ch = chAll + (R_xlen_t)t * nEdge;
      graph_geodesic_phylo(&nTip, &nNode, par, ch, &nEdge, &allNodes,
                           interim.data());
      const double inv = 1.0 / (t + 1);
      R_xlen_t p = 0;
      for (int j = 0; j < nTip - 1; ++j) {
        const int colOffset = allNodes * j;
        for (int i = j + 1; i < nTip; ++i, ++p) {
          const double x = lg[interim[i + colOffset]];
          const double delta = x - mean[p];
          mean[p] += delta * inv;
          m2[p] += delta * (x - mean[p]);
        }
      }
    }
    const double denom = nTree > 1 ? nTree - 1 : 1;
    for (R_xlen_t p = 0; p < nPairs; ++p) {
      aves[p] = mean[p];
      const double var = m2[p] / denom;
      const double sd = var > 0 ? std::sqrt(var) : 0.0;
      devs[p] = std::isfinite(sd) ? sd : 0.0;
    }
  } else {
    /* median and/or mad: buffer distances (tree-major, contiguous writes),
     * then reduce each pair with nth_element. */
    std::vector<double> buf(nPairs * (R_xlen_t)nTree);
    for (int t = 0; t < nTree; ++t) {
      R_CheckUserInterrupt();
      const int *par = parAll + (R_xlen_t)t * nEdge;
      const int *ch = chAll + (R_xlen_t)t * nEdge;
      graph_geodesic_phylo(&nTip, &nNode, par, ch, &nEdge, &allNodes,
                           interim.data());
      double *col = &buf[(R_xlen_t)t * nPairs];
      R_xlen_t p = 0;
      for (int j = 0; j < nTip - 1; ++j) {
        const int colOffset = allNodes * j;
        for (int i = j + 1; i < nTip; ++i, ++p) {
          col[p] = lg[interim[i + colOffset]];
        }
      }
    }

    // Scratch buffers reused across pairs (no per-pair heap churn).
    std::vector<double> row(nTree), scratch(nTree), absDev(nTree);
    const double denom = nTree > 1 ? nTree - 1 : 1;
    const bool needSum = (whichAve == 1 || whichDev == 1);  // mean or sd
    const bool needMed = (whichAve == 2 || whichDev == 2);  // median or mad
    for (R_xlen_t p = 0; p < nPairs; ++p) {
      for (int t = 0; t < nTree; ++t) {
        row[t] = buf[(R_xlen_t)t * nPairs + p];
      }

      double mean = 0.0, med = 0.0;
      if (needSum) {
        double sum = 0.0;
        for (int t = 0; t < nTree; ++t) sum += row[t];
        mean = sum / nTree;
      }
      if (needMed) {
        // The MAD centre is the same row median used by average = "median",
        // so compute it once and share it.
        std::copy(row.begin(), row.end(), scratch.begin());
        med = median_inplace(scratch.data(), scratch.data() + nTree);
      }

      aves[p] = (whichAve == 1) ? mean : med;

      double dev;
      if (whichDev == 1) {
        double ss = 0.0;
        for (int t = 0; t < nTree; ++t) {
          const double d = row[t] - mean;
          ss += d * d;
        }
        const double var = ss / denom;
        dev = var > 0 ? std::sqrt(var) : 0.0;
      } else {
        for (int t = 0; t < nTree; ++t) absDev[t] = std::fabs(row[t] - med);
        dev = 1.4826 * median_inplace(absDev.data(), absDev.data() + nTree);
      }
      devs[p] = std::isfinite(dev) ? dev : 0.0;
    }
  }

  /* meanAve = mean deviation-average over all pairs (normalising constant). */
  double sumAve = 0.0;
  for (R_xlen_t p = 0; p < nPairs; ++p) sumAve += aves[p];
  const double meanAve = sumAve / nPairs;

  /* Fold each pair's deviation onto both its leaves -- equivalent to
   * rowMeans of the symmetric (zero-diagonal) deviation matrix, without
   * building it. */
  SEXP result = PROTECT(allocVector(REALSXP, nTip));
  double *res = REAL(result);
  for (int i = 0; i < nTip; ++i) res[i] = 0.0;
  {
    R_xlen_t p = 0;
    for (int j = 0; j < nTip - 1; ++j) {
      for (int i = j + 1; i < nTip; ++i, ++p) {
        res[i] += devs[p];
        res[j] += devs[p];
      }
    }
  }
  const double norm = 1.0 / ((double)nTip * meanAve);
  for (int i = 0; i < nTip; ++i) res[i] *= norm;

  UNPROTECT(1);
  return result;
}
