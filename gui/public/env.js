window.nomadEnv = {
  "appBase": "https://nomad-lab.eu/prod/v1",
  "northBase": "http://localhost:9000/fairdi/nomad/latest/north",
  "keycloakBase": "https://nomad-lab.eu/fairdi/keycloak/auth/",
  "keycloakRealm": "fairdi_nomad_prod",
  "keycloakClientId": "nomad_public",
  "debug": false,
  "encyclopediaBase": "https://nomad-lab.eu/prod/rae/encyclopedia/#",
  "aitoolkitEnabled": false,
  "oasis": false,
  "version": {},
  "globalLoginRequired": false,
  "servicesUploadLimit": 10,
  "ui": {
    "theme": {
      "title": "NOMAD"
    },
    "unit_systems": {
      "selected": "Custom"
    },
    "entry": {
      "cards": {
        "exclude": [
          "relatedResources"
        ],
        "options": {
          "sections": {
            "error": "Could not render section card."
          },
          "definitions": {
            "error": "Could not render definitions card."
          },
          "nexus": {
            "error": "Could not render NeXus card."
          },
          "material": {
            "error": "Could not render material card."
          },
          "solarcell": {
            "error": "Could not render solar cell properties."
          },
          "electronic": {
            "error": "Could not render electronic properties."
          },
          "vibrational": {
            "error": "Could not render vibrational properties."
          },
          "mechanical": {
            "error": "Could not render mechanical properties."
          },
          "thermodynamic": {
            "error": "Could not render thermodynamic properties."
          },
          "structural": {
            "error": "Could not render structural properties."
          },
          "dynamical": {
            "error": "Could not render dynamical properties."
          },
          "geometry_optimization": {
            "error": "Could not render geometry optimization."
          },
          "eels": {
            "error": "Could not render EELS properties."
          },
          "workflow": {
            "error": "Could not render workflow card."
          },
          "references": {
            "error": "Could not render references card."
          },
          "relatedResources": {
            "error": "Could not render related resources card."
          }
        }
      }
    },
    "apps": {
      "options": {
        "entries": {
          "label": "Entries",
          "path": "entries",
          "resource": "entries",
          "breadcrumb": "Entries",
          "category": "All",
          "description": "Search entries across all domains",
          "help": {
            "title": "Entries search",
            "content": "This page allows you to search **entries** within NOMAD.\nEntries represent any individual data items that have\nbeen uploaded to NOMAD, no matter whether they come from\ntheoretical calculations, experiments, lab notebooks or\nany other source of data. This allows you to perform\ncross-domain queries, but if you are interested in a\nspecific subfield, you should see if a specific\napplication exists for it in the explore menu to get\nmore details."
          },
          "pagination": {
            "order_by": "upload_create_time",
            "order": "desc",
            "page_size": 20
          },
          "columns": {
            "options": {
              "entry_name": {
                "label": "Name",
                "align": "left"
              },
              "results.material.chemical_formula_hill": {
                "label": "Formula",
                "align": "left"
              },
              "entry_type": {
                "label": "Entry type",
                "align": "left"
              },
              "upload_create_time": {
                "label": "Upload time",
                "align": "left"
              },
              "authors": {
                "label": "Authors",
                "align": "left"
              },
              "results.method.method_name": {
                "label": "Method name",
                "align": "left"
              },
              "results.method.simulation.program_name": {
                "label": "Program name",
                "align": "left"
              },
              "results.method.simulation.dft.basis_set_name": {
                "label": "Basis set name",
                "align": "left"
              },
              "results.method.simulation.dft.xc_functional_type": {
                "label": "XC Functional Type",
                "align": "left"
              },
              "results.material.structural_type": {
                "label": "Dimensionality",
                "align": "left"
              },
              "results.material.symmetry.crystal_system": {
                "label": "Crystal system",
                "align": "left"
              },
              "results.material.symmetry.space_group_symbol": {
                "label": "Space group symbol",
                "align": "left"
              },
              "results.material.symmetry.space_group_number": {
                "label": "Space group number",
                "align": "left"
              },
              "results.eln.lab_ids": {
                "label": "Lab IDs",
                "align": "left"
              },
              "results.eln.sections": {
                "label": "Sections",
                "align": "left"
              },
              "results.eln.methods": {
                "label": "Methods",
                "align": "left"
              },
              "results.eln.tags": {
                "label": "Tags",
                "align": "left"
              },
              "results.eln.instruments": {
                "label": "Instruments",
                "align": "left"
              },
              "mainfile": {
                "label": "Mainfile",
                "align": "left"
              },
              "comment": {
                "label": "Comment",
                "align": "left"
              },
              "references": {
                "label": "References",
                "align": "left"
              },
              "datasets": {
                "label": "Datasets",
                "align": "left"
              },
              "published": {
                "label": "Access",
                "align": "left"
              }
            },
            "selected": [
              "entry_name",
              "results.material.chemical_formula_hill",
              "entry_type",
              "upload_create_time",
              "authors"
            ]
          },
          "rows": {
            "actions": {
              "enabled": true
            },
            "details": {
              "enabled": true
            },
            "selection": {
              "enabled": true
            }
          },
          "filter_menus": {
            "options": {
              "material": {
                "label": "Material",
                "level": 0,
                "size": "small"
              },
              "elements": {
                "label": "Elements / Formula",
                "level": 1,
                "size": "large"
              },
              "structure": {
                "label": "Structure",
                "level": 1,
                "size": "small"
              },
              "method": {
                "label": "Method",
                "level": 0,
                "size": "small"
              },
              "dft": {
                "label": "DFT",
                "level": 1,
                "size": "small"
              },
              "gw": {
                "label": "GW",
                "level": 1,
                "size": "small"
              },
              "projection": {
                "label": "Projection",
                "level": 1,
                "size": "small"
              },
              "dmft": {
                "label": "DMFT",
                "level": 1,
                "size": "small"
              },
              "eels": {
                "label": "EELS",
                "level": 1,
                "size": "small"
              },
              "workflow": {
                "label": "Workflow",
                "level": 0,
                "size": "small"
              },
              "molecular_dynamics": {
                "label": "Molecular dynamics",
                "level": 1,
                "size": "small"
              },
              "geometry_optimization": {
                "label": "Geometry Optimization",
                "level": 1,
                "size": "small"
              },
              "properties": {
                "label": "Properties",
                "level": 0,
                "size": "small"
              },
              "electronic": {
                "label": "Electronic",
                "level": 1,
                "size": "small"
              },
              "vibrational": {
                "label": "Vibrational",
                "level": 1,
                "size": "small"
              },
              "mechanical": {
                "label": "Mechanical",
                "level": 1,
                "size": "small"
              },
              "usecases": {
                "label": "Use Cases",
                "level": 0,
                "size": "small"
              },
              "solarcell": {
                "label": "Solar Cells",
                "level": 1,
                "size": "small"
              },
              "author": {
                "label": "Author / Origin / Dataset",
                "level": 0,
                "size": "medium"
              },
              "metadata": {
                "label": "Visibility / IDs / Schema",
                "level": 0,
                "size": "small"
              },
              "optimade": {
                "label": "Optimade",
                "level": 0,
                "size": "medium"
              }
            }
          },
          "filters": {
            "exclude": [
              "mainfile",
              "entry_name",
              "combine"
            ]
          }
        },
        "calculations": {
          "label": "Calculations",
          "path": "calculations",
          "resource": "entries",
          "breadcrumb": "Calculations",
          "category": "Theory",
          "description": "Search calculations",
          "help": {
            "title": "Calculations",
            "content": "This page allows you to search **calculations** within\nNOMAD.  Calculations typically come from a specific\nsimulation software that uses an approximate model to\ninvestigate and report different physical properties."
          },
          "pagination": {
            "order_by": "upload_create_time",
            "order": "desc",
            "page_size": 20
          },
          "columns": {
            "options": {
              "results.material.chemical_formula_hill": {
                "label": "Formula",
                "align": "left"
              },
              "results.method.simulation.program_name": {
                "label": "Program name",
                "align": "left"
              },
              "results.method.method_name": {
                "label": "Method name",
                "align": "left"
              },
              "upload_create_time": {
                "label": "Upload time",
                "align": "left"
              },
              "authors": {
                "label": "Authors",
                "align": "left"
              },
              "results.method.simulation.dft.basis_set_name": {
                "label": "Basis set name",
                "align": "left"
              },
              "results.method.simulation.dft.xc_functional_type": {
                "label": "XC Functional Type",
                "align": "left"
              },
              "results.material.structural_type": {
                "label": "Dimensionality",
                "align": "left"
              },
              "results.material.symmetry.crystal_system": {
                "label": "Crystal system",
                "align": "left"
              },
              "results.material.symmetry.space_group_symbol": {
                "label": "Space group symbol",
                "align": "left"
              },
              "results.material.symmetry.space_group_number": {
                "label": "Space group number",
                "align": "left"
              },
              "results.eln.lab_ids": {
                "label": "Lab IDs",
                "align": "left"
              },
              "results.eln.sections": {
                "label": "Sections",
                "align": "left"
              },
              "results.eln.methods": {
                "label": "Methods",
                "align": "left"
              },
              "results.eln.tags": {
                "label": "Tags",
                "align": "left"
              },
              "results.eln.instruments": {
                "label": "Instruments",
                "align": "left"
              },
              "entry_name": {
                "label": "Name",
                "align": "left"
              },
              "entry_type": {
                "label": "Entry type",
                "align": "left"
              },
              "mainfile": {
                "label": "Mainfile",
                "align": "left"
              },
              "comment": {
                "label": "Comment",
                "align": "left"
              },
              "references": {
                "label": "References",
                "align": "left"
              },
              "datasets": {
                "label": "Datasets",
                "align": "left"
              },
              "published": {
                "label": "Access",
                "align": "left"
              }
            },
            "selected": [
              "results.material.chemical_formula_hill",
              "results.method.simulation.program_name",
              "results.method.method_name",
              "upload_create_time",
              "authors"
            ]
          },
          "rows": {
            "actions": {
              "enabled": true
            },
            "details": {
              "enabled": true
            },
            "selection": {
              "enabled": true
            }
          },
          "filter_menus": {
            "options": {
              "material": {
                "label": "Material",
                "level": 0,
                "size": "small"
              },
              "elements": {
                "label": "Elements / Formula",
                "level": 1,
                "size": "large"
              },
              "structure": {
                "label": "Structure",
                "level": 1,
                "size": "small"
              },
              "method": {
                "label": "Method",
                "level": 0,
                "size": "small"
              },
              "dft": {
                "label": "DFT",
                "level": 1,
                "size": "small"
              },
              "gw": {
                "label": "GW",
                "level": 1,
                "size": "small"
              },
              "projection": {
                "label": "Projection",
                "level": 1,
                "size": "small"
              },
              "dmft": {
                "label": "DMFT",
                "level": 1,
                "size": "small"
              },
              "workflow": {
                "label": "Workflow",
                "level": 0,
                "size": "small"
              },
              "molecular_dynamics": {
                "label": "Molecular dynamics",
                "level": 1,
                "size": "small"
              },
              "geometry_optimization": {
                "label": "Geometry Optimization",
                "level": 1,
                "size": "small"
              },
              "properties": {
                "label": "Properties",
                "level": 0,
                "size": "small"
              },
              "electronic": {
                "label": "Electronic",
                "level": 1,
                "size": "small"
              },
              "vibrational": {
                "label": "Vibrational",
                "level": 1,
                "size": "small"
              },
              "mechanical": {
                "label": "Mechanical",
                "level": 1,
                "size": "small"
              },
              "author": {
                "label": "Author / Origin / Dataset",
                "level": 0,
                "size": "medium"
              },
              "metadata": {
                "label": "Visibility / IDs / Schema",
                "level": 0,
                "size": "small"
              },
              "optimade": {
                "label": "Optimade",
                "level": 0,
                "size": "medium"
              }
            }
          },
          "filters": {
            "exclude": [
              "mainfile",
              "entry_name",
              "combine"
            ]
          },
          "filters_locked": {
            "quantities": "results.method.simulation.program_name"
          }
        },
        "materials": {
          "label": "Materials",
          "path": "materials",
          "resource": "materials",
          "breadcrumb": "Materials",
          "category": "Theory",
          "description": "Search materials that are identified from calculations",
          "help": {
            "title": "Materials",
            "content": "This page allows you to search **materials** within\nNOMAD. NOMAD can often automatically detect the material\nfrom individual calculations that contain the full\natomistic structure and can then group the data by using\nthese detected materials. This allows you to search\nindividual materials which have properties that are\naggregated from several entries. Following the link for\na specific material will take you to the corresponding\n[NOMAD Encyclopedia](https://nomad-lab.eu/prod/rae/encyclopedia/#/search)\npage for that material. NOMAD Encyclopedia is a service\nthat is specifically oriented towards materials property\nexploration.\n\nNotice that by default the properties that you search\ncan be combined from several different entries. If\ninstead you wish to search for a material with an\nindividual entry fullfilling your search criteria,\nuncheck the **combine results from several\nentries**-checkbox."
          },
          "pagination": {
            "order_by": "chemical_formula_hill",
            "order": "asc",
            "page_size": 20
          },
          "columns": {
            "options": {
              "chemical_formula_hill": {
                "label": "Formula",
                "align": "left"
              },
              "structural_type": {
                "label": "Dimensionality",
                "align": "left"
              },
              "symmetry.structure_name": {
                "label": "Structure name",
                "align": "left"
              },
              "symmetry.space_group_number": {
                "label": "Space group number",
                "align": "left"
              },
              "symmetry.crystal_system": {
                "label": "Crystal system",
                "align": "left"
              },
              "symmetry.space_group_symbol": {
                "label": "Space group symbol",
                "align": "left"
              },
              "material_id": {
                "label": "Material ID",
                "align": "left"
              }
            },
            "selected": [
              "chemical_formula_hill",
              "structural_type",
              "symmetry.structure_name",
              "symmetry.space_group_number",
              "symmetry.crystal_system"
            ]
          },
          "rows": {
            "actions": {
              "enabled": true
            },
            "details": {
              "enabled": false
            },
            "selection": {
              "enabled": false
            }
          },
          "filter_menus": {
            "options": {
              "material": {
                "label": "Material",
                "level": 0,
                "size": "small"
              },
              "elements": {
                "label": "Elements / Formula",
                "level": 1,
                "size": "large"
              },
              "structure": {
                "label": "Structure",
                "level": 1,
                "size": "small"
              },
              "method": {
                "label": "Method",
                "level": 0,
                "size": "small"
              },
              "dft": {
                "label": "DFT",
                "level": 1,
                "size": "small"
              },
              "gw": {
                "label": "GW",
                "level": 1,
                "size": "small"
              },
              "projection": {
                "label": "Projection",
                "level": 1,
                "size": "small"
              },
              "dmft": {
                "label": "DMFT",
                "level": 1,
                "size": "small"
              },
              "workflow": {
                "label": "Workflow",
                "level": 0,
                "size": "small"
              },
              "molecular_dynamics": {
                "label": "Molecular dynamics",
                "level": 1,
                "size": "small"
              },
              "geometry_optimization": {
                "label": "Geometry Optimization",
                "level": 1,
                "size": "small"
              },
              "properties": {
                "label": "Properties",
                "level": 0,
                "size": "small"
              },
              "electronic": {
                "label": "Electronic",
                "level": 1,
                "size": "small"
              },
              "vibrational": {
                "label": "Vibrational",
                "level": 1,
                "size": "small"
              },
              "mechanical": {
                "label": "Mechanical",
                "level": 1,
                "size": "small"
              },
              "author": {
                "label": "Author / Origin / Dataset",
                "level": 0,
                "size": "medium"
              },
              "metadata": {
                "label": "Visibility / IDs / Schema",
                "level": 0,
                "size": "small"
              },
              "optimade": {
                "label": "Optimade",
                "level": 0,
                "size": "medium"
              },
              "combine": {
                "level": 0,
                "size": "small",
                "actions": {
                  "options": {
                    "combine": {
                      "type": "checkbox",
                      "label": "Combine results from several entries",
                      "quantity": "combine"
                    }
                  }
                }
              }
            }
          },
          "filters": {
            "exclude": [
              "mainfile",
              "entry_name"
            ]
          }
        },
        "eln": {
          "label": "ELN",
          "path": "eln",
          "resource": "entries",
          "breadcrumb": "ELN",
          "category": "Experiment",
          "description": "Search electronic lab notebooks",
          "help": {
            "title": "ELN search",
            "content": "This page allows you to specifically seach **electronic\nlab notebooks (ELNs)** within NOMAD.  It is very similar\nto the entries search, but with a reduced filter set and\nspecialized arrangement of default columns."
          },
          "pagination": {
            "order_by": "upload_create_time",
            "order": "desc",
            "page_size": 20
          },
          "columns": {
            "options": {
              "entry_name": {
                "label": "Name",
                "align": "left"
              },
              "entry_type": {
                "label": "Entry type",
                "align": "left"
              },
              "upload_create_time": {
                "label": "Upload time",
                "align": "left"
              },
              "authors": {
                "label": "Authors",
                "align": "left"
              },
              "results.material.chemical_formula_hill": {
                "label": "Formula",
                "align": "left"
              },
              "results.method.method_name": {
                "label": "Method name",
                "align": "left"
              },
              "results.eln.lab_ids": {
                "label": "Lab IDs",
                "align": "left"
              },
              "results.eln.sections": {
                "label": "Sections",
                "align": "left"
              },
              "results.eln.methods": {
                "label": "Methods",
                "align": "left"
              },
              "results.eln.tags": {
                "label": "Tags",
                "align": "left"
              },
              "results.eln.instruments": {
                "label": "Instruments",
                "align": "left"
              },
              "mainfile": {
                "label": "Mainfile",
                "align": "left"
              },
              "comment": {
                "label": "Comment",
                "align": "left"
              },
              "references": {
                "label": "References",
                "align": "left"
              },
              "datasets": {
                "label": "Datasets",
                "align": "left"
              },
              "published": {
                "label": "Access",
                "align": "left"
              }
            },
            "selected": [
              "entry_name",
              "entry_type",
              "upload_create_time",
              "authors"
            ]
          },
          "rows": {
            "actions": {
              "enabled": true
            },
            "details": {
              "enabled": true
            },
            "selection": {
              "enabled": true
            }
          },
          "filter_menus": {
            "options": {
              "material": {
                "label": "Material",
                "level": 0,
                "size": "small"
              },
              "elements": {
                "label": "Elements / Formula",
                "level": 1,
                "size": "large"
              },
              "eln": {
                "label": "Electronic Lab Notebook",
                "level": 0,
                "size": "small"
              },
              "custom_quantities": {
                "label": "User Defined Quantities",
                "level": 0,
                "size": "large"
              },
              "author": {
                "label": "Author / Origin / Dataset",
                "level": 0,
                "size": "medium"
              },
              "metadata": {
                "label": "Visibility / IDs / Schema",
                "level": 0,
                "size": "small"
              },
              "optimade": {
                "label": "Optimade",
                "level": 0,
                "size": "medium"
              }
            }
          },
          "filters": {
            "exclude": [
              "mainfile",
              "entry_name",
              "combine"
            ]
          },
          "filters_locked": {
            "quantities": "data"
          }
        },
        "eels": {
          "label": "EELS",
          "path": "eels",
          "resource": "entries",
          "breadcrumb": "EELS",
          "category": "Experiment",
          "description": "Search electron energy loss spectroscopy experiments",
          "help": {
            "title": "EELS",
            "content": "This page allows you to spefically search **Electron\nEnergy Loss Spectroscopy (EELS) experiments** within\nNOMAD. It is similar to the entries search, but with a\nreduced filter set and specialized arrangement of\ndefault columns."
          },
          "pagination": {
            "order_by": "upload_create_time",
            "order": "desc",
            "page_size": 20
          },
          "columns": {
            "options": {
              "results.material.chemical_formula_hill": {
                "label": "Formula",
                "align": "left"
              },
              "results.properties.spectroscopy.eels.detector_type": {
                "label": "Detector type",
                "align": "left"
              },
              "results.properties.spectroscopy.eels.resolution": {
                "label": "Resolution",
                "align": "left"
              },
              "upload_create_time": {
                "label": "Upload time",
                "align": "left"
              },
              "authors": {
                "label": "Authors",
                "align": "left"
              },
              "results.properties.spectroscopy.eels.min_energy": {
                "align": "left"
              },
              "results.properties.spectroscopy.eels.max_energy": {
                "align": "left"
              },
              "entry_name": {
                "label": "Name",
                "align": "left"
              },
              "entry_type": {
                "label": "Entry type",
                "align": "left"
              },
              "mainfile": {
                "label": "Mainfile",
                "align": "left"
              },
              "comment": {
                "label": "Comment",
                "align": "left"
              },
              "references": {
                "label": "References",
                "align": "left"
              },
              "datasets": {
                "label": "Datasets",
                "align": "left"
              },
              "published": {
                "label": "Access",
                "align": "left"
              }
            },
            "selected": [
              "results.material.chemical_formula_hill",
              "results.properties.spectroscopy.eels.detector_type",
              "results.properties.spectroscopy.eels.resolution",
              "upload_create_time",
              "authors"
            ]
          },
          "rows": {
            "actions": {
              "enabled": true
            },
            "details": {
              "enabled": true
            },
            "selection": {
              "enabled": true
            }
          },
          "filter_menus": {
            "options": {
              "material": {
                "label": "Material",
                "level": 0,
                "size": "small"
              },
              "elements": {
                "label": "Elements / Formula",
                "level": 1,
                "size": "large"
              },
              "method": {
                "label": "Method",
                "level": 0,
                "size": "small"
              },
              "eels": {
                "label": "EELS",
                "level": 1,
                "size": "small"
              },
              "author": {
                "label": "Author / Origin / Dataset",
                "level": 0,
                "size": "medium"
              },
              "metadata": {
                "label": "Visibility / IDs / Schema",
                "level": 0,
                "size": "small"
              },
              "optimade": {
                "label": "Optimade",
                "level": 0,
                "size": "medium"
              }
            }
          },
          "filters": {
            "exclude": [
              "mainfile",
              "entry_name",
              "combine"
            ]
          },
          "filters_locked": {
            "results.method.method_name": "EELS"
          }
        },
        "solarcells": {
          "label": "Solar Cells",
          "path": "solarcells",
          "resource": "entries",
          "breadcrumb": "Solar Cells",
          "category": "Use Cases",
          "description": "Search solar cells",
          "help": {
            "title": "Solar cells",
            "content": "This page allows you to search **solar cell data**\nwithin NOMAD. The filter menu on the left and the shown\ndefault columns are specifically designed for solar cell\nexploration. The dashboard directly shows useful\ninteractive statistics about the data."
          },
          "pagination": {
            "order_by": "results.properties.optoelectronic.solar_cell.efficiency",
            "order": "desc",
            "page_size": 20
          },
          "columns": {
            "options": {
              "results.material.chemical_formula_descriptive": {
                "label": "Descriptive Formula",
                "align": "left"
              },
              "results.properties.optoelectronic.solar_cell.efficiency": {
                "label": "Efficiency (%)",
                "align": "left",
                "format": {
                  "decimals": 2,
                  "mode": "standard"
                }
              },
              "results.properties.optoelectronic.solar_cell.open_circuit_voltage": {
                "label": "Open circuit voltage",
                "align": "left",
                "unit": "V",
                "format": {
                  "decimals": 3,
                  "mode": "standard"
                }
              },
              "results.properties.optoelectronic.solar_cell.short_circuit_current_density": {
                "label": "Short circuit current density",
                "align": "left",
                "unit": "A/m**2",
                "format": {
                  "decimals": 3,
                  "mode": "standard"
                }
              },
              "results.properties.optoelectronic.solar_cell.fill_factor": {
                "label": "Fill factor",
                "align": "left",
                "format": {
                  "decimals": 3,
                  "mode": "standard"
                }
              },
              "references": {
                "label": "References",
                "align": "left"
              },
              "results.material.chemical_formula_hill": {
                "label": "Formula",
                "align": "left"
              },
              "results.material.structural_type": {
                "label": "Dimensionality",
                "align": "left"
              },
              "results.properties.optoelectronic.solar_cell.illumination_intensity": {
                "label": "Illum. intensity",
                "align": "left",
                "unit": "W/m**2",
                "format": {
                  "decimals": 3,
                  "mode": "standard"
                }
              },
              "results.eln.lab_ids": {
                "label": "Lab IDs",
                "align": "left"
              },
              "results.eln.sections": {
                "label": "Sections",
                "align": "left"
              },
              "results.eln.methods": {
                "label": "Methods",
                "align": "left"
              },
              "results.eln.tags": {
                "label": "Tags",
                "align": "left"
              },
              "results.eln.instruments": {
                "label": "Instruments",
                "align": "left"
              },
              "entry_name": {
                "label": "Name",
                "align": "left"
              },
              "entry_type": {
                "label": "Entry type",
                "align": "left"
              },
              "mainfile": {
                "label": "Mainfile",
                "align": "left"
              },
              "upload_create_time": {
                "label": "Upload time",
                "align": "left"
              },
              "authors": {
                "label": "Authors",
                "align": "left"
              },
              "comment": {
                "label": "Comment",
                "align": "left"
              },
              "datasets": {
                "label": "Datasets",
                "align": "left"
              },
              "published": {
                "label": "Access",
                "align": "left"
              }
            },
            "selected": [
              "results.material.chemical_formula_descriptive",
              "results.properties.optoelectronic.solar_cell.efficiency",
              "results.properties.optoelectronic.solar_cell.open_circuit_voltage",
              "results.properties.optoelectronic.solar_cell.short_circuit_current_density",
              "results.properties.optoelectronic.solar_cell.fill_factor",
              "references"
            ]
          },
          "rows": {
            "actions": {
              "enabled": true
            },
            "details": {
              "enabled": true
            },
            "selection": {
              "enabled": true
            }
          },
          "filter_menus": {
            "options": {
              "material": {
                "label": "Absorber Material",
                "level": 0,
                "size": "small"
              },
              "elements": {
                "label": "Elements / Formula",
                "level": 1,
                "size": "large"
              },
              "structure": {
                "label": "Structure",
                "level": 1,
                "size": "small"
              },
              "electronic": {
                "label": "Electronic Properties",
                "level": 0,
                "size": "small"
              },
              "solarcell": {
                "label": "Solar Cell Properties",
                "level": 0,
                "size": "small"
              },
              "eln": {
                "label": "Electronic Lab Notebook",
                "level": 0,
                "size": "small"
              },
              "custom_quantities": {
                "label": "User Defined Quantities",
                "level": 0,
                "size": "large"
              },
              "author": {
                "label": "Author / Origin / Dataset",
                "level": 0,
                "size": "medium"
              },
              "metadata": {
                "label": "Visibility / IDs / Schema",
                "level": 0,
                "size": "small"
              },
              "optimade": {
                "label": "Optimade",
                "level": 0,
                "size": "medium"
              }
            }
          },
          "filters": {
            "exclude": [
              "mainfile",
              "entry_name",
              "combine"
            ]
          },
          "dashboard": {
            "widgets": [
              {
                "type": "periodictable",
                "layout": {
                  "xxl": {
                    "minH": 8,
                    "minW": 12,
                    "h": 8,
                    "w": 13,
                    "x": 0,
                    "y": 0
                  },
                  "xl": {
                    "minH": 8,
                    "minW": 12,
                    "h": 8,
                    "w": 12,
                    "x": 0,
                    "y": 0
                  },
                  "lg": {
                    "minH": 8,
                    "minW": 12,
                    "h": 8,
                    "w": 12,
                    "x": 0,
                    "y": 0
                  },
                  "md": {
                    "minH": 8,
                    "minW": 12,
                    "h": 8,
                    "w": 12,
                    "x": 0,
                    "y": 0
                  },
                  "sm": {
                    "minH": 8,
                    "minW": 12,
                    "h": 8,
                    "w": 12,
                    "x": 0,
                    "y": 16
                  }
                },
                "quantity": "results.material.elements",
                "scale": "linear"
              },
              {
                "type": "scatterplot",
                "layout": {
                  "xxl": {
                    "minH": 3,
                    "minW": 3,
                    "h": 8,
                    "w": 12,
                    "x": 24,
                    "y": 0
                  },
                  "xl": {
                    "minH": 3,
                    "minW": 3,
                    "h": 8,
                    "w": 9,
                    "x": 12,
                    "y": 0
                  },
                  "lg": {
                    "minH": 3,
                    "minW": 3,
                    "h": 6,
                    "w": 12,
                    "x": 0,
                    "y": 8
                  },
                  "md": {
                    "minH": 3,
                    "minW": 3,
                    "h": 6,
                    "w": 9,
                    "x": 0,
                    "y": 8
                  },
                  "sm": {
                    "minH": 3,
                    "minW": 3,
                    "h": 5,
                    "w": 6,
                    "x": 0,
                    "y": 0
                  }
                },
                "x": "results.properties.optoelectronic.solar_cell.open_circuit_voltage",
                "y": "results.properties.optoelectronic.solar_cell.efficiency",
                "color": "results.properties.optoelectronic.solar_cell.short_circuit_current_density",
                "size": 1000,
                "autorange": true
              },
              {
                "type": "scatterplot",
                "layout": {
                  "xxl": {
                    "minH": 3,
                    "minW": 3,
                    "h": 8,
                    "w": 11,
                    "x": 13,
                    "y": 0
                  },
                  "xl": {
                    "minH": 3,
                    "minW": 3,
                    "h": 8,
                    "w": 9,
                    "x": 21,
                    "y": 0
                  },
                  "lg": {
                    "minH": 3,
                    "minW": 3,
                    "h": 6,
                    "w": 12,
                    "x": 0,
                    "y": 14
                  },
                  "md": {
                    "minH": 3,
                    "minW": 3,
                    "h": 6,
                    "w": 9,
                    "x": 9,
                    "y": 8
                  },
                  "sm": {
                    "minH": 3,
                    "minW": 3,
                    "h": 5,
                    "w": 6,
                    "x": 6,
                    "y": 0
                  }
                },
                "x": "results.properties.optoelectronic.solar_cell.open_circuit_voltage",
                "y": "results.properties.optoelectronic.solar_cell.efficiency",
                "color": "results.properties.optoelectronic.solar_cell.device_architecture",
                "size": 1000,
                "autorange": true
              },
              {
                "type": "terms",
                "layout": {
                  "xxl": {
                    "minH": 3,
                    "minW": 3,
                    "h": 6,
                    "w": 6,
                    "x": 14,
                    "y": 8
                  },
                  "xl": {
                    "minH": 3,
                    "minW": 3,
                    "h": 6,
                    "w": 6,
                    "x": 14,
                    "y": 8
                  },
                  "lg": {
                    "minH": 3,
                    "minW": 3,
                    "h": 6,
                    "w": 6,
                    "x": 12,
                    "y": 0
                  },
                  "md": {
                    "minH": 3,
                    "minW": 3,
                    "h": 4,
                    "w": 6,
                    "x": 12,
                    "y": 4
                  },
                  "sm": {
                    "minH": 3,
                    "minW": 3,
                    "h": 6,
                    "w": 4,
                    "x": 0,
                    "y": 10
                  }
                },
                "quantity": "results.properties.optoelectronic.solar_cell.device_stack",
                "scale": "linear",
                "showinput": true
              },
              {
                "type": "histogram",
                "layout": {
                  "xxl": {
                    "minH": 3,
                    "minW": 3,
                    "h": 3,
                    "w": 8,
                    "x": 0,
                    "y": 8
                  },
                  "xl": {
                    "minH": 3,
                    "minW": 3,
                    "h": 3,
                    "w": 8,
                    "x": 0,
                    "y": 11
                  },
                  "lg": {
                    "minH": 3,
                    "minW": 3,
                    "h": 4,
                    "w": 12,
                    "x": 12,
                    "y": 12
                  },
                  "md": {
                    "minH": 3,
                    "minW": 3,
                    "h": 3,
                    "w": 8,
                    "x": 10,
                    "y": 17
                  },
                  "sm": {
                    "minH": 3,
                    "minW": 3,
                    "h": 3,
                    "w": 8,
                    "x": 4,
                    "y": 13
                  }
                },
                "quantity": "results.properties.optoelectronic.solar_cell.illumination_intensity",
                "scale": "1/4",
                "autorange": true,
                "showinput": true,
                "nbins": 30
              },
              {
                "type": "terms",
                "layout": {
                  "xxl": {
                    "minH": 3,
                    "minW": 3,
                    "h": 6,
                    "w": 6,
                    "x": 8,
                    "y": 8
                  },
                  "xl": {
                    "minH": 3,
                    "minW": 3,
                    "h": 6,
                    "w": 6,
                    "x": 8,
                    "y": 8
                  },
                  "lg": {
                    "minH": 3,
                    "minW": 3,
                    "h": 6,
                    "w": 6,
                    "x": 18,
                    "y": 0
                  },
                  "md": {
                    "minH": 3,
                    "minW": 3,
                    "h": 4,
                    "w": 6,
                    "x": 12,
                    "y": 0
                  },
                  "sm": {
                    "minH": 3,
                    "minW": 3,
                    "h": 5,
                    "w": 4,
                    "x": 0,
                    "y": 5
                  }
                },
                "quantity": "results.properties.optoelectronic.solar_cell.absorber_fabrication",
                "scale": "linear",
                "showinput": true
              },
              {
                "type": "histogram",
                "layout": {
                  "xxl": {
                    "minH": 3,
                    "minW": 8,
                    "h": 3,
                    "w": 8,
                    "x": 0,
                    "y": 11
                  },
                  "xl": {
                    "minH": 3,
                    "minW": 8,
                    "h": 3,
                    "w": 8,
                    "x": 0,
                    "y": 8
                  },
                  "lg": {
                    "minH": 3,
                    "minW": 8,
                    "h": 4,
                    "w": 12,
                    "x": 12,
                    "y": 16
                  },
                  "md": {
                    "minH": 3,
                    "minW": 8,
                    "h": 3,
                    "w": 8,
                    "x": 10,
                    "y": 14
                  },
                  "sm": {
                    "minH": 3,
                    "minW": 8,
                    "h": 3,
                    "w": 8,
                    "x": 4,
                    "y": 10
                  }
                },
                "quantity": "results.properties.electronic.band_structure_electronic.band_gap.value",
                "scale": "1/4",
                "autorange": false,
                "showinput": false,
                "nbins": 30
              },
              {
                "type": "terms",
                "layout": {
                  "xxl": {
                    "minH": 3,
                    "minW": 3,
                    "h": 6,
                    "w": 6,
                    "x": 20,
                    "y": 8
                  },
                  "xl": {
                    "minH": 3,
                    "minW": 3,
                    "h": 6,
                    "w": 5,
                    "x": 25,
                    "y": 8
                  },
                  "lg": {
                    "minH": 3,
                    "minW": 3,
                    "h": 6,
                    "w": 6,
                    "x": 18,
                    "y": 6
                  },
                  "md": {
                    "minH": 3,
                    "minW": 3,
                    "h": 6,
                    "w": 5,
                    "x": 0,
                    "y": 14
                  },
                  "sm": {
                    "minH": 3,
                    "minW": 3,
                    "h": 5,
                    "w": 4,
                    "x": 4,
                    "y": 5
                  }
                },
                "quantity": "results.properties.optoelectronic.solar_cell.electron_transport_layer",
                "scale": "linear",
                "showinput": true
              },
              {
                "type": "terms",
                "layout": {
                  "xxl": {
                    "minH": 3,
                    "minW": 3,
                    "h": 6,
                    "w": 6,
                    "x": 26,
                    "y": 8
                  },
                  "xl": {
                    "minH": 3,
                    "minW": 3,
                    "h": 6,
                    "w": 5,
                    "x": 20,
                    "y": 8
                  },
                  "lg": {
                    "minH": 3,
                    "minW": 3,
                    "h": 6,
                    "w": 6,
                    "x": 12,
                    "y": 6
                  },
                  "md": {
                    "minH": 3,
                    "minW": 3,
                    "h": 6,
                    "w": 5,
                    "x": 5,
                    "y": 14
                  },
                  "sm": {
                    "minH": 3,
                    "minW": 3,
                    "h": 5,
                    "w": 4,
                    "x": 8,
                    "y": 5
                  }
                },
                "quantity": "results.properties.optoelectronic.solar_cell.hole_transport_layer",
                "scale": "linear",
                "showinput": true
              }
            ]
          },
          "filters_locked": {
            "sections": "nomad.datamodel.results.SolarCell"
          }
        }
      }
    },
    "north": {},
    "example_uploads": {}
  }
}
