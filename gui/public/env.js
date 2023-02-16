window.nomadEnv = {
  "appBase": "http://localhost:8000/fairdi/nomad/latest",
  "northBase": "http://localhost:9000/fairdi/nomad/latest/north",
  "keycloakBase": "https://nomad-lab.eu/fairdi/keycloak/auth/",
  "keycloakRealm": "fairdi_nomad_test",
  "keycloakClientId": "nomad_public",
  "debug": false,
  "encyclopediaBase": "https://nomad-lab.eu/prod/rae/encyclopedia/#",
  "aitoolkitEnabled": false,
  "oasis": false,
  "version": {},
  "globalLoginRequired": false,
  "servicesUploadLimit": 10,
  "ui": {
    "entry_context": {
      "overview": {
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
    "search_contexts": {
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
            "enable": [
              "entry_name",
              "results.material.chemical_formula_hill",
              "entry_type",
              "upload_create_time",
              "authors"
            ],
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
                "label": "Method name"
              },
              "results.method.simulation.program_name": {
                "label": "Program name"
              },
              "results.method.simulation.dft.basis_set_name": {
                "label": "Basis set name"
              },
              "results.method.simulation.dft.xc_functional_type": {
                "label": "XC Functional Type"
              },
              "results.material.structural_type": {
                "label": "Dimensionality"
              },
              "results.material.symmetry.crystal_system": {
                "label": "Crystal system"
              },
              "results.material.symmetry.space_group_symbol": {
                "label": "Space group symbol"
              },
              "results.material.symmetry.space_group_number": {
                "label": "Space group number"
              },
              "results.eln.lab_ids": {
                "label": "Lab IDs"
              },
              "results.eln.sections": {
                "label": "Sections"
              },
              "results.eln.methods": {
                "label": "Methods"
              },
              "results.eln.tags": {
                "label": "Tags"
              },
              "results.eln.instruments": {
                "label": "Instruments"
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
                "label": "Access"
              }
            }
          },
          "rows": {
            "actions": {
              "enable": true
            },
            "details": {
              "enable": true
            },
            "selection": {
              "enable": true
            }
          },
          "filter_menus": {
            "options": {
              "material": {
                "label": "Material",
                "level": 0
              },
              "elements": {
                "label": "Elements / Formula",
                "level": 1,
                "size": "large",
                "menu_items": {}
              },
              "structure": {
                "label": "Structure",
                "level": 1,
                "size": "small",
                "menu_items": {}
              },
              "method": {
                "label": "Method",
                "level": 0,
                "menu_items": {}
              },
              "dft": {
                "label": "DFT",
                "level": 1,
                "size": "small",
                "menu_items": {}
              },
              "gw": {
                "label": "GW",
                "level": 1,
                "size": "small",
                "menu_items": {}
              },
              "projection": {
                "label": "Projection",
                "level": 1,
                "size": "small",
                "menu_items": {}
              },
              "dmft": {
                "label": "DMFT",
                "level": 1,
                "size": "small",
                "menu_items": {}
              },
              "eels": {
                "label": "EELS",
                "level": 1,
                "size": "small",
                "menu_items": {}
              },
              "workflow": {
                "label": "Workflow",
                "level": 0
              },
              "molecular_dynamics": {
                "label": "Molecular dynamics",
                "level": 1,
                "size": "small",
                "menu_items": {}
              },
              "geometry_optimization": {
                "label": "Geometry Optimization",
                "level": 1,
                "size": "small",
                "menu_items": {}
              },
              "properties": {
                "label": "Properties",
                "level": 0
              },
              "electronic": {
                "label": "Electronic",
                "level": 1,
                "size": "small",
                "menu_items": {}
              },
              "vibrational": {
                "label": "Vibrational",
                "level": 1,
                "size": "small",
                "menu_items": {}
              },
              "mechanical": {
                "label": "Mechanical",
                "level": 1,
                "size": "small",
                "menu_items": {}
              },
              "usecases": {
                "label": "Use Cases",
                "level": 0,
                "size": "small"
              },
              "solarcell": {
                "label": "Solar Cells",
                "level": 1,
                "size": "small",
                "menu_items": {}
              },
              "author": {
                "label": "Author / Origin / Dataset",
                "level": 0,
                "size": "medium",
                "menu_items": {}
              },
              "metadata": {
                "label": "Visibility / IDs / Schema",
                "level": 0,
                "size": "small",
                "menu_items": {}
              },
              "optimade": {
                "label": "Optimade",
                "level": 0,
                "size": "medium",
                "menu_items": {}
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
            "enable": [
              "results.material.chemical_formula_hill",
              "results.method.simulation.program_name",
              "results.method.method_name",
              "upload_create_time",
              "authors"
            ],
            "options": {
              "results.material.chemical_formula_hill": {
                "label": "Formula",
                "align": "left"
              },
              "results.method.simulation.program_name": {
                "label": "Program name"
              },
              "results.method.method_name": {
                "label": "Method name"
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
                "label": "Basis set name"
              },
              "results.method.simulation.dft.xc_functional_type": {
                "label": "XC Functional Type"
              },
              "results.material.structural_type": {
                "label": "Dimensionality"
              },
              "results.material.symmetry.crystal_system": {
                "label": "Crystal system"
              },
              "results.material.symmetry.space_group_symbol": {
                "label": "Space group symbol"
              },
              "results.material.symmetry.space_group_number": {
                "label": "Space group number"
              },
              "results.eln.lab_ids": {
                "label": "Lab IDs"
              },
              "results.eln.sections": {
                "label": "Sections"
              },
              "results.eln.methods": {
                "label": "Methods"
              },
              "results.eln.tags": {
                "label": "Tags"
              },
              "results.eln.instruments": {
                "label": "Instruments"
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
                "label": "Access"
              }
            }
          },
          "rows": {
            "actions": {
              "enable": true
            },
            "details": {
              "enable": true
            },
            "selection": {
              "enable": true
            }
          },
          "filter_menus": {
            "options": {
              "material": {
                "label": "Material",
                "level": 0
              },
              "elements": {
                "label": "Elements / Formula",
                "level": 1,
                "size": "large",
                "menu_items": {}
              },
              "structure": {
                "label": "Structure",
                "level": 1,
                "size": "small",
                "menu_items": {}
              },
              "method": {
                "label": "Method",
                "level": 0,
                "menu_items": {}
              },
              "dft": {
                "label": "DFT",
                "level": 1,
                "size": "small",
                "menu_items": {}
              },
              "gw": {
                "label": "GW",
                "level": 1,
                "size": "small",
                "menu_items": {}
              },
              "projection": {
                "label": "Projection",
                "level": 1,
                "size": "small",
                "menu_items": {}
              },
              "dmft": {
                "label": "DMFT",
                "level": 1,
                "size": "small",
                "menu_items": {}
              },
              "workflow": {
                "label": "Workflow",
                "level": 0
              },
              "molecular_dynamics": {
                "label": "Molecular dynamics",
                "level": 1,
                "size": "small",
                "menu_items": {}
              },
              "geometry_optimization": {
                "label": "Geometry Optimization",
                "level": 1,
                "size": "small",
                "menu_items": {}
              },
              "properties": {
                "label": "Properties",
                "level": 0
              },
              "electronic": {
                "label": "Electronic",
                "level": 1,
                "size": "small",
                "menu_items": {}
              },
              "vibrational": {
                "label": "Vibrational",
                "level": 1,
                "size": "small",
                "menu_items": {}
              },
              "mechanical": {
                "label": "Mechanical",
                "level": 1,
                "size": "small",
                "menu_items": {}
              },
              "author": {
                "label": "Author / Origin / Dataset",
                "level": 0,
                "size": "medium",
                "menu_items": {}
              },
              "metadata": {
                "label": "Visibility / IDs / Schema",
                "level": 0,
                "size": "small",
                "menu_items": {}
              },
              "optimade": {
                "label": "Optimade",
                "level": 0,
                "size": "medium",
                "menu_items": {}
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
            "order": "asc"
          },
          "columns": {
            "enable": [
              "chemical_formula_hill",
              "structural_type",
              "symmetry.structure_name",
              "symmetry.space_group_number",
              "symmetry.crystal_system"
            ],
            "options": {
              "chemical_formula_hill": {
                "label": "Formula",
                "align": "left"
              },
              "structural_type": {
                "label": "Dimensionality"
              },
              "symmetry.structure_name": {
                "label": "Structure name"
              },
              "symmetry.space_group_number": {
                "label": "Space group number"
              },
              "symmetry.crystal_system": {
                "label": "Crystal system"
              },
              "symmetry.space_group_symbol": {
                "label": "Space group symbol"
              },
              "material_id": {
                "label": "Material ID"
              }
            }
          },
          "rows": {
            "actions": {
              "enable": true
            },
            "details": {
              "enable": false
            },
            "selection": {
              "enable": false
            }
          },
          "filter_menus": {
            "options": {
              "material": {
                "label": "Material",
                "level": 0
              },
              "elements": {
                "label": "Elements / Formula",
                "level": 1,
                "size": "large",
                "menu_items": {}
              },
              "structure": {
                "label": "Structure",
                "level": 1,
                "size": "small",
                "menu_items": {}
              },
              "method": {
                "label": "Method",
                "level": 0,
                "menu_items": {}
              },
              "dft": {
                "label": "DFT",
                "level": 1,
                "size": "small",
                "menu_items": {}
              },
              "gw": {
                "label": "GW",
                "level": 1,
                "size": "small",
                "menu_items": {}
              },
              "projection": {
                "label": "Projection",
                "level": 1,
                "size": "small",
                "menu_items": {}
              },
              "dmft": {
                "label": "DMFT",
                "level": 1,
                "size": "small",
                "menu_items": {}
              },
              "workflow": {
                "label": "Workflow",
                "level": 0
              },
              "molecular_dynamics": {
                "label": "Molecular dynamics",
                "level": 1,
                "size": "small",
                "menu_items": {}
              },
              "geometry_optimization": {
                "label": "Geometry Optimization",
                "level": 1,
                "size": "small",
                "menu_items": {}
              },
              "properties": {
                "label": "Properties",
                "level": 0,
                "size": "small"
              },
              "electronic": {
                "label": "Electronic",
                "level": 1,
                "size": "small",
                "menu_items": {}
              },
              "vibrational": {
                "label": "Vibrational",
                "level": 1,
                "size": "small",
                "menu_items": {}
              },
              "mechanical": {
                "label": "Mechanical",
                "level": 1,
                "size": "small",
                "menu_items": {}
              },
              "author": {
                "label": "Author / Origin / Dataset",
                "level": 0,
                "size": "medium",
                "menu_items": {}
              },
              "metadata": {
                "label": "Visibility / IDs / Schema",
                "level": 0,
                "size": "small",
                "menu_items": {}
              },
              "optimade": {
                "label": "Optimade",
                "level": 0,
                "size": "medium",
                "menu_items": {}
              },
              "combine": {
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
            "enable": [
              "entry_name",
              "entry_type",
              "upload_create_time",
              "authors"
            ],
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
                "label": "Method name"
              },
              "results.eln.lab_ids": {
                "label": "Lab IDs"
              },
              "results.eln.sections": {
                "label": "Sections"
              },
              "results.eln.methods": {
                "label": "Methods"
              },
              "results.eln.tags": {
                "label": "Tags"
              },
              "results.eln.instruments": {
                "label": "Instruments"
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
                "label": "Access"
              }
            }
          },
          "rows": {
            "actions": {
              "enable": true
            },
            "details": {
              "enable": true
            },
            "selection": {
              "enable": true
            }
          },
          "filter_menus": {
            "options": {
              "material": {
                "label": "Material",
                "level": 0
              },
              "elements": {
                "label": "Elements / Formula",
                "level": 1,
                "size": "large",
                "menu_items": {}
              },
              "eln": {
                "label": "Electronic Lab Notebook",
                "level": 0,
                "size": "small",
                "menu_items": {}
              },
              "custom_quantities": {
                "label": "User Defined Quantities",
                "level": 0,
                "size": "large",
                "menu_items": {}
              },
              "author": {
                "label": "Author / Origin / Dataset",
                "level": 0,
                "size": "medium",
                "menu_items": {}
              },
              "metadata": {
                "label": "Visibility / IDs / Schema",
                "level": 0,
                "size": "small",
                "menu_items": {}
              },
              "optimade": {
                "label": "Optimade",
                "level": 0,
                "size": "medium",
                "menu_items": {}
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
            "enable": [
              "results.material.chemical_formula_hill",
              "results.properties.spectroscopy.eels.detector_type",
              "results.properties.spectroscopy.eels.resolution",
              "upload_create_time",
              "authors"
            ],
            "options": {
              "results.material.chemical_formula_hill": {
                "label": "Formula",
                "align": "left"
              },
              "results.properties.spectroscopy.eels.detector_type": {
                "label": "Detector type"
              },
              "results.properties.spectroscopy.eels.resolution": {
                "label": "Resolution"
              },
              "upload_create_time": {
                "label": "Upload time",
                "align": "left"
              },
              "authors": {
                "label": "Authors",
                "align": "left"
              },
              "results.properties.spectroscopy.eels.min_energy": {},
              "results.properties.spectroscopy.eels.max_energy": {},
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
                "label": "Access"
              }
            }
          },
          "rows": {
            "actions": {
              "enable": true
            },
            "details": {
              "enable": true
            },
            "selection": {
              "enable": true
            }
          },
          "filter_menus": {
            "options": {
              "material": {
                "label": "Material",
                "level": 0
              },
              "elements": {
                "label": "Elements / Formula",
                "level": 1,
                "size": "large",
                "menu_items": {}
              },
              "method": {
                "label": "Method",
                "level": 0,
                "size": "small"
              },
              "eels": {
                "label": "EELS",
                "level": 1,
                "size": "small",
                "menu_items": {}
              },
              "author": {
                "label": "Author / Origin / Dataset",
                "level": 0,
                "size": "medium",
                "menu_items": {}
              },
              "metadata": {
                "label": "Visibility / IDs / Schema",
                "level": 0,
                "size": "small",
                "menu_items": {}
              },
              "optimade": {
                "label": "Optimade",
                "level": 0,
                "size": "medium",
                "menu_items": {}
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
          "dashboard": {
            "widgets": [
              {
                "type": "periodictable",
                "scale": "linear",
                "quantity": "results.material.elements",
                "layout": {
                  "xxl": {
                    "minH": 8,
                    "minW": 12,
                    "h": 8,
                    "w": 13,
                    "y": 0,
                    "x": 0
                  },
                  "xl": {
                    "minH": 8,
                    "minW": 12,
                    "h": 8,
                    "w": 12,
                    "y": 0,
                    "x": 0
                  },
                  "lg": {
                    "minH": 8,
                    "minW": 12,
                    "h": 8,
                    "w": 12,
                    "y": 0,
                    "x": 0
                  },
                  "md": {
                    "minH": 8,
                    "minW": 12,
                    "h": 8,
                    "w": 12,
                    "y": 0,
                    "x": 0
                  },
                  "sm": {
                    "minH": 8,
                    "minW": 12,
                    "h": 8,
                    "w": 12,
                    "y": 16,
                    "x": 0
                  }
                }
              },
              {
                "type": "scatterplot",
                "autorange": true,
                "size": 1000,
                "color": "results.properties.optoelectronic.solar_cell.short_circuit_current_density",
                "y": "results.properties.optoelectronic.solar_cell.efficiency",
                "x": "results.properties.optoelectronic.solar_cell.open_circuit_voltage",
                "layout": {
                  "xxl": {
                    "minH": 3,
                    "minW": 3,
                    "h": 8,
                    "w": 12,
                    "y": 0,
                    "x": 24
                  },
                  "xl": {
                    "minH": 3,
                    "minW": 3,
                    "h": 8,
                    "w": 9,
                    "y": 0,
                    "x": 12
                  },
                  "lg": {
                    "minH": 3,
                    "minW": 3,
                    "h": 6,
                    "w": 12,
                    "y": 8,
                    "x": 0
                  },
                  "md": {
                    "minH": 3,
                    "minW": 3,
                    "h": 6,
                    "w": 9,
                    "y": 8,
                    "x": 0
                  },
                  "sm": {
                    "minH": 3,
                    "minW": 3,
                    "h": 5,
                    "w": 6,
                    "y": 0,
                    "x": 0
                  }
                }
              },
              {
                "type": "scatterplot",
                "autorange": true,
                "size": 1000,
                "color": "results.properties.optoelectronic.solar_cell.device_architecture",
                "y": "results.properties.optoelectronic.solar_cell.efficiency",
                "x": "results.properties.optoelectronic.solar_cell.open_circuit_voltage",
                "layout": {
                  "xxl": {
                    "minH": 3,
                    "minW": 3,
                    "h": 8,
                    "w": 11,
                    "y": 0,
                    "x": 13
                  },
                  "xl": {
                    "minH": 3,
                    "minW": 3,
                    "h": 8,
                    "w": 9,
                    "y": 0,
                    "x": 21
                  },
                  "lg": {
                    "minH": 3,
                    "minW": 3,
                    "h": 6,
                    "w": 12,
                    "y": 14,
                    "x": 0
                  },
                  "md": {
                    "minH": 3,
                    "minW": 3,
                    "h": 6,
                    "w": 9,
                    "y": 8,
                    "x": 9
                  },
                  "sm": {
                    "minH": 3,
                    "minW": 3,
                    "h": 5,
                    "w": 6,
                    "y": 0,
                    "x": 6
                  }
                }
              },
              {
                "type": "terms",
                "inputfields": true,
                "scale": "linear",
                "quantity": "results.properties.optoelectronic.solar_cell.device_stack",
                "layout": {
                  "xxl": {
                    "minH": 3,
                    "minW": 3,
                    "h": 6,
                    "w": 6,
                    "y": 8,
                    "x": 14
                  },
                  "xl": {
                    "minH": 3,
                    "minW": 3,
                    "h": 6,
                    "w": 6,
                    "y": 8,
                    "x": 14
                  },
                  "lg": {
                    "minH": 3,
                    "minW": 3,
                    "h": 6,
                    "w": 6,
                    "y": 0,
                    "x": 12
                  },
                  "md": {
                    "minH": 3,
                    "minW": 3,
                    "h": 4,
                    "w": 6,
                    "y": 4,
                    "x": 12
                  },
                  "sm": {
                    "minH": 3,
                    "minW": 3,
                    "h": 6,
                    "w": 4,
                    "y": 10,
                    "x": 0
                  }
                }
              },
              {
                "type": "histogram",
                "autorange": true,
                "nbins": 30,
                "scale": "1/4",
                "quantity": "results.properties.optoelectronic.solar_cell.illumination_intensity",
                "layout": {
                  "xxl": {
                    "minH": 3,
                    "minW": 3,
                    "h": 3,
                    "w": 8,
                    "y": 8,
                    "x": 0
                  },
                  "xl": {
                    "minH": 3,
                    "minW": 3,
                    "h": 3,
                    "w": 8,
                    "y": 11,
                    "x": 0
                  },
                  "lg": {
                    "minH": 3,
                    "minW": 3,
                    "h": 4,
                    "w": 12,
                    "y": 12,
                    "x": 12
                  },
                  "md": {
                    "minH": 3,
                    "minW": 3,
                    "h": 3,
                    "w": 8,
                    "y": 17,
                    "x": 10
                  },
                  "sm": {
                    "minH": 3,
                    "minW": 3,
                    "h": 3,
                    "w": 8,
                    "y": 13,
                    "x": 4
                  }
                }
              },
              {
                "type": "terms",
                "inputfields": true,
                "scale": "linear",
                "quantity": "results.properties.optoelectronic.solar_cell.absorber_fabrication",
                "layout": {
                  "xxl": {
                    "minH": 3,
                    "minW": 3,
                    "h": 6,
                    "w": 6,
                    "y": 8,
                    "x": 8
                  },
                  "xl": {
                    "minH": 3,
                    "minW": 3,
                    "h": 6,
                    "w": 6,
                    "y": 8,
                    "x": 8
                  },
                  "lg": {
                    "minH": 3,
                    "minW": 3,
                    "h": 6,
                    "w": 6,
                    "y": 0,
                    "x": 18
                  },
                  "md": {
                    "minH": 3,
                    "minW": 3,
                    "h": 4,
                    "w": 6,
                    "y": 0,
                    "x": 12
                  },
                  "sm": {
                    "minH": 3,
                    "minW": 3,
                    "h": 5,
                    "w": 4,
                    "y": 5,
                    "x": 0
                  }
                }
              },
              {
                "type": "histogram",
                "inputfields": false,
                "autorange": false,
                "nbins": 30,
                "scale": "1/4",
                "quantity": "results.properties.electronic.band_structure_electronic.band_gap.value",
                "layout": {
                  "xxl": {
                    "minH": 3,
                    "minW": 8,
                    "h": 3,
                    "w": 8,
                    "y": 11,
                    "x": 0
                  },
                  "xl": {
                    "minH": 3,
                    "minW": 8,
                    "h": 3,
                    "w": 8,
                    "y": 8,
                    "x": 0
                  },
                  "lg": {
                    "minH": 3,
                    "minW": 8,
                    "h": 4,
                    "w": 12,
                    "y": 16,
                    "x": 12
                  },
                  "md": {
                    "minH": 3,
                    "minW": 8,
                    "h": 3,
                    "w": 8,
                    "y": 14,
                    "x": 10
                  },
                  "sm": {
                    "minH": 3,
                    "minW": 8,
                    "h": 3,
                    "w": 8,
                    "y": 10,
                    "x": 4
                  }
                }
              },
              {
                "type": "terms",
                "inputfields": true,
                "scale": "linear",
                "quantity": "results.properties.optoelectronic.solar_cell.electron_transport_layer",
                "layout": {
                  "xxl": {
                    "minH": 3,
                    "minW": 3,
                    "h": 6,
                    "w": 6,
                    "y": 8,
                    "x": 20
                  },
                  "xl": {
                    "minH": 3,
                    "minW": 3,
                    "h": 6,
                    "w": 5,
                    "y": 8,
                    "x": 25
                  },
                  "lg": {
                    "minH": 3,
                    "minW": 3,
                    "h": 6,
                    "w": 6,
                    "y": 6,
                    "x": 18
                  },
                  "md": {
                    "minH": 3,
                    "minW": 3,
                    "h": 6,
                    "w": 5,
                    "y": 14,
                    "x": 0
                  },
                  "sm": {
                    "minH": 3,
                    "minW": 3,
                    "h": 5,
                    "w": 4,
                    "y": 5,
                    "x": 4
                  }
                }
              },
              {
                "type": "terms",
                "inputfields": true,
                "scale": "linear",
                "quantity": "results.properties.optoelectronic.solar_cell.hole_transport_layer",
                "layout": {
                  "xxl": {
                    "minH": 3,
                    "minW": 3,
                    "h": 6,
                    "w": 6,
                    "y": 8,
                    "x": 26
                  },
                  "xl": {
                    "minH": 3,
                    "minW": 3,
                    "h": 6,
                    "w": 5,
                    "y": 8,
                    "x": 20
                  },
                  "lg": {
                    "minH": 3,
                    "minW": 3,
                    "h": 6,
                    "w": 6,
                    "y": 6,
                    "x": 12
                  },
                  "md": {
                    "minH": 3,
                    "minW": 3,
                    "h": 6,
                    "w": 5,
                    "y": 14,
                    "x": 5
                  },
                  "sm": {
                    "minH": 3,
                    "minW": 3,
                    "h": 5,
                    "w": 4,
                    "y": 5,
                    "x": 8
                  }
                }
              }
            ]
          },
          "columns": {
            "enable": [
              "results.material.chemical_formula_descriptive",
              "results.properties.optoelectronic.solar_cell.efficiency",
              "results.properties.optoelectronic.solar_cell.open_circuit_voltage",
              "results.properties.optoelectronic.solar_cell.short_circuit_current_density",
              "results.properties.optoelectronic.solar_cell.fill_factor",
              "references"
            ],
            "options": {
              "results.material.chemical_formula_descriptive": {
                "label": "Descriptive Formula",
                "align": "left"
              },
              "results.properties.optoelectronic.solar_cell.efficiency": {
                "label": "Efficiency (%)",
                "format": {
                  "decimals": 2,
                  "mode": "standard"
                }
              },
              "results.properties.optoelectronic.solar_cell.open_circuit_voltage": {
                "label": "Open circuit voltage",
                "unit": "V",
                "format": {
                  "decimals": 3,
                  "mode": "standard"
                }
              },
              "results.properties.optoelectronic.solar_cell.short_circuit_current_density": {
                "label": "Short circuit current density",
                "unit": "A/m**2",
                "format": {
                  "decimals": 3,
                  "mode": "standard"
                }
              },
              "results.properties.optoelectronic.solar_cell.fill_factor": {
                "label": "Fill factor",
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
                "label": "Dimensionality"
              },
              "results.properties.optoelectronic.solar_cell.illumination_intensity": {
                "label": "Illum. intensity",
                "unit": "W/m**2",
                "format": {
                  "decimals": 3,
                  "mode": "standard"
                }
              },
              "results.eln.lab_ids": {
                "label": "Lab IDs"
              },
              "results.eln.sections": {
                "label": "Sections"
              },
              "results.eln.methods": {
                "label": "Methods"
              },
              "results.eln.tags": {
                "label": "Tags"
              },
              "results.eln.instruments": {
                "label": "Instruments"
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
                "label": "Access"
              }
            }
          },
          "rows": {
            "actions": {
              "enable": true
            },
            "details": {
              "enable": true
            },
            "selection": {
              "enable": true
            }
          },
          "filter_menus": {
            "options": {
              "material": {
                "label": "Absorber Material",
                "level": 0
              },
              "elements": {
                "label": "Elements / Formula",
                "level": 1,
                "size": "large",
                "menu_items": {}
              },
              "structure": {
                "label": "Structure",
                "level": 1,
                "size": "small",
                "menu_items": {}
              },
              "electronic": {
                "label": "Electronic Properties",
                "level": 0,
                "size": "small",
                "menu_items": {}
              },
              "solarcell": {
                "label": "Solar Cell Properties",
                "level": 0,
                "size": "small",
                "menu_items": {}
              },
              "eln": {
                "label": "Electronic Lab Notebook",
                "level": 0,
                "size": "small",
                "menu_items": {}
              },
              "custom_quantities": {
                "label": "User Defined Quantities",
                "level": 0,
                "size": "large",
                "menu_items": {}
              },
              "author": {
                "label": "Author / Origin / Dataset",
                "level": 0,
                "size": "medium",
                "menu_items": {}
              },
              "metadata": {
                "label": "Visibility / IDs / Schema",
                "level": 0,
                "size": "small",
                "menu_items": {}
              },
              "optimade": {
                "label": "Optimade",
                "level": 0,
                "size": "medium",
                "menu_items": {}
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
            "sections": "nomad.datamodel.results.SolarCell"
          }
        }
      }
    },
    "example_uploads": {
      "include": null,
      "exclude": null
    },
    "default_unit_system": "Custom",
    "north_enabled": false,
    "theme": {
      "title": "NOMAD"
    }
  }
}
