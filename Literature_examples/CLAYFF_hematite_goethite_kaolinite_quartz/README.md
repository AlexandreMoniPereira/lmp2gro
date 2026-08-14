# Molecular Models of Hematite, Goethite, Kaolinite, and Quartz: Surface Terminations, Ionic Interactions, Nano-topography, and Water Coordination

This guide provides a walkthrough on how to use **lmp2gro** to convert LAMMPS data files for hematite, goethite, kaolinite, and quartz into GROMACS format. These models were first discussed and characterized in the following publication:

[![Journal Article](https://img.shields.io/badge/DOI-10.1016/j.colsurfa.2022.129585-blue)](https://doi.org/10.1016/j.colsurfa.2022.129585)

The complete dataset is also available free of charge on Zenodo:

[![Zenodo Dataset](https://img.shields.io/badge/DOI-10.5281/zenodo.6685648-blue)](https://doi.org/10.5281/zenodo.6685648)

---

## LAMMPS Data Files

The data files are provided herein. Detailed information regarding the construction, force field parameters, and structural validation for each mineral model can be found in the original paper and the Zenodo repository.

## Conversion Process

To convert any data file in this folder using `lmp2gro`, use the following command syntax:

```bash
lmp2gro <data_file> -r <RESNAME> --folder <output_folder>
```

For example, the aluminum-terminated kaolinite (`KAO_Al`) was generated using:
```bash
lmp2gro kao-Al-0.data -r KAO --folder KAO_Al
```

## Citations

If you utilize these files in your research, please cite the **lmp2gro** tool, the **Clayff** force field, and the specific mineral models. Below are the BibTeX entries for your convenience:


```bibtex
@article{Pereira2026_lmp2gro,
    author = {Moni Pereira, Alexandre and da Silva Martins, Jarede and Albuquerque, Marcelo and Costa, Luciano T.},
    title = {lmp2gro: A Python Tool for Converting LAMMPS Data Files into GROMACS Topologies},
    journal = {Journal of Chemical Information and Modeling},
    volume = {66},
    number = {15},
    pages = {8698-8704},
    year = {2026},
    month = {07},
    issn = {1549-9596},
    doi = {10.1021/acs.jcim.6c01430},
    url = {https://doi.org/10.1021/acs.jcim.6c01430},
    eprint = {https://pubs.acs.org/jcisd8/article-pdf/66/15/8698/66054108/acs.jcim.6c01430.pdf},
}

@article{Cygan2004_Clayff,
    title = {Molecular Models of Hydroxide, Oxyhydroxide, and Clay Phases and the Development of a General Force Field},
    volume = {108},
    number = {4},
    journal = {The Journal of Physical Chemistry B},
    doi = {10.1021/jp0363287},
    author = {Cygan, Randall T. and Liang, Jian-Jie and Kalinichev, Andrey G.},
    year = {2004},
    pages = {1255–1266}
}

@article{Filippov2022_Molecular_Models,
    title = {Molecular models of hematite, goethite, kaolinite, and quartz: Surface terminations, ionic interactions, nano topography, and water coordination},
    volume = {650},
    journal = {Colloids and Surfaces A: Physicochemical and Engineering Aspects},
    doi = {10.1016/j.colsurfa.2022.129585},
    author = {Filippov, Lev O. and Silva, Lucas A. and Pereira, Alexandre M. and Bastos, Leonardo C. and Correia, Julio C.G. and Silva, Klaydison and Pi\c{c}arra, Alexandre and Foucaud, Yann},
    year = {2022},
    pages = {129585}
}

@misc{Filippov2022_zenodo_models,
    doi = {10.5281/ZENODO.6685649},
    author = {Silva, Lucas A. and Pereira, Alexandre M. and Correia, Julio C. G. and Bastos, Leonardo C. and Silva, Klaydison and Pi\c{c}arra, Alexandre and Foucaud, Yann and Filippov, Lev O.},
    title = {Molecular Models of hematite, goethite, kaolinite, and quartz surfaces},
    publisher = {Zenodo},
    year = {2022}
}
```

