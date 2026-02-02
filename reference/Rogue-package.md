# Rogue

"Rogue" implements approaches to identify rogue taxa in phylogenetic
analysis. Rogues are wildcard leaves whose uncertain position, perhaps a
result of missing or conflicting data, reduces the resolution of
consensus trees (Kearney 2002) . Consensus trees that omit rogue taxa
can be more informative.

## Details

"Rogue" allows the user to select a concept of "information" by which
the quality of consensus trees should be evaluated, and a heuristic
approach by which rogue taxa should be identified.

Rogue detection using the phylogenetic and clustering information
content measures (SPIC, SCIC) (Smith 2022) is implemented using a quick
heuristic that drops the least "stable" leaves one at a time, using an
*ad hoc* definition of stability (Smith 2022) ; and by a more exhaustive
(and time-consuming) approach that considers dropping all possible sets
of up to *n* leaves (Aberer et al. 2013) .

The latter heuristic is implemented for the relative bipartition
"information" content and Pattengale's criterion *via*
[RogueNaRok](https://cme.h-its.org/exelixis/web/software/roguenarok/roguenarok.html)
(Aberer et al. 2013) .

### Citing "Rogue"

If you find this package useful in your work, Please consider citing
Smith (2021).

To cite the underlying methods, please cite Aberer et al. (2013)
("RogueNaRok") or Smith (2022) (SPIC, SCIC), as appropriate.

## References

Aberer AJ, Krompass D, Stamatakis A (2013). “Pruning rogue taxa improves
phylogenetic accuracy: an efficient algorithm and webservice.”
*Systematic Biology*, **62**(1), 162–166.
[doi:10.1093/sysbio/sys078](https://doi.org/10.1093/sysbio/sys078) .  
  
Kearney M (2002). “Fragmentary taxa, missing data, and ambiguity:
mistaken assumptions and conclusions.” *Systematic Biology*, **51**(2),
369–381.
[doi:10.1080/10635150252899824](https://doi.org/10.1080/10635150252899824)
.  
  
Smith MR (2022). “Using information theory to detect rogue taxa and
improve consensus trees.” *Systematic Biology*, **71**(5), 986–1008.
[doi:10.1093/sysbio/syab099](https://doi.org/10.1093/sysbio/syab099) .

## See also

Useful links:

- <https://github.com/ms609/Rogue/>

- <https://github.com/aberer/RogueNaRok/>

- <https://github.com/ms609/RogueNaRok/>

- Report bugs at <https://github.com/ms609/Rogue/issues/>

## Author

**Maintainer**: Martin R. Smith <martin.smith@durham.ac.uk>
([ORCID](https://orcid.org/0000-0001-5660-1727)) \[copyright holder\]

Authors:

- Andre J. Aberer <andre.aberer@googlemail.com> (RogueNaRok) \[copyright
  holder\]
