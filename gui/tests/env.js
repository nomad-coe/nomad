window.nomadEnv = {
  "appBase": "http://localhost:8000/fairdi/nomad/latest",
  "northBase": "http://localhost:9000/fairdi/nomad/latest/north",
  "keycloakBase": "https://nomad-lab.eu/fairdi/keycloak/auth/",
  "keycloakRealm": "fairdi_nomad_test",
  "keycloakClientId": "nomad_public",
  "debug": false,
  "encyclopediaBase": "https://nomad-lab.eu/prod/rae/encyclopedia/#",
  "aitoolkitEnabled": true,
  "oasis": false,
  "version": {},
  "globalLoginRequired": false,
  "servicesUploadLimit": 10,
  "appTokenMaxExpiresIn": 2592000,
  "ui": {
    "app_base": "http://localhost:8000/fairdi/nomad/latest",
    "north_base": "http://localhost:9000/fairdi/nomad/latest/north",
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
          "heterogeneouscatalyst": {
            "error": "Could not render catalyst properties."
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
          "spectroscopic": {
            "error": "Could not render spectroscopic properties."
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
      "exclude": [
        "heterogeneouscatalyst"
      ],
      "options": {
        "entries": {
          "label": "Entries",
          "path": "entries",
          "resource": "entries",
          "category": "All",
          "description": "Search entries across all domains",
          "readme": "This page allows you to search **entries** within NOMAD.\nEntries represent any individual data items that have\nbeen uploaded to NOMAD, no matter whether they come from\ntheoretical calculations, experiments, lab notebooks or\nany other source of data. This allows you to perform\ncross-domain queries, but if you are interested in a\nspecific subfield, you should see if a specific\napplication exists for it in the explore menu to get\nmore details.",
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
              "upload_name": {
                "label": "Upload name",
                "align": "left"
              },
              "upload_id": {
                "label": "Upload id",
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
              "results.method.simulation.dft.xc_functional_type": {
                "label": "XC Functional Type",
                "align": "left"
              },
              "results.method.simulation.precision.apw_cutoff": {
                "label": "APW Cutoff",
                "align": "left"
              },
              "results.method.simulation.precision.basis_set": {
                "label": "Basis Set",
                "align": "left"
              },
              "results.method.simulation.precision.k_line_density": {
                "label": "k-line Density",
                "align": "left"
              },
              "results.method.simulation.precision.native_tier": {
                "label": "Code-specific tier",
                "align": "left"
              },
              "results.method.simulation.precision.planewave_cutoff": {
                "label": "Plane-wave cutoff",
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
                "size": "s"
              },
              "elements": {
                "label": "Elements / Formula",
                "level": 1,
                "size": "xl"
              },
              "structure": {
                "label": "Structure / Symmetry",
                "level": 1,
                "size": "s"
              },
              "method": {
                "label": "Method",
                "level": 0,
                "size": "s"
              },
              "precision": {
                "label": "Precision",
                "level": 1,
                "size": "s"
              },
              "dft": {
                "label": "DFT",
                "level": 1,
                "size": "s"
              },
              "gw": {
                "label": "GW",
                "level": 1,
                "size": "s"
              },
              "bse": {
                "label": "BSE",
                "level": 1,
                "size": "s"
              },
              "projection": {
                "label": "Projection",
                "level": 1,
                "size": "s"
              },
              "dmft": {
                "label": "DMFT",
                "level": 1,
                "size": "s"
              },
              "eels": {
                "label": "EELS",
                "level": 1,
                "size": "s"
              },
              "workflow": {
                "label": "Workflow",
                "level": 0,
                "size": "s"
              },
              "molecular_dynamics": {
                "label": "Molecular dynamics",
                "level": 1,
                "size": "s"
              },
              "geometry_optimization": {
                "label": "Geometry Optimization",
                "level": 1,
                "size": "s"
              },
              "properties": {
                "label": "Properties",
                "level": 0,
                "size": "s"
              },
              "electronic": {
                "label": "Electronic",
                "level": 1,
                "size": "s"
              },
              "vibrational": {
                "label": "Vibrational",
                "level": 1,
                "size": "s"
              },
              "mechanical": {
                "label": "Mechanical",
                "level": 1,
                "size": "s"
              },
              "usecases": {
                "label": "Use Cases",
                "level": 0,
                "size": "s"
              },
              "solarcell": {
                "label": "Solar Cells",
                "level": 1,
                "size": "s"
              },
              "heterogeneouscatalyst": {
                "label": "Heterogeneous Catalysis",
                "level": 1,
                "size": "s"
              },
              "author": {
                "label": "Author / Origin / Dataset",
                "level": 0,
                "size": "m"
              },
              "metadata": {
                "label": "Visibility / IDs / Schema",
                "level": 0,
                "size": "s"
              },
              "optimade": {
                "label": "Optimade",
                "level": 0,
                "size": "m"
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
          "category": "Theory",
          "description": "Search calculations",
          "readme": "This page allows you to search **calculations** within\nNOMAD. Calculations typically come from a specific\nsimulation software that uses an approximate model to\ninvestigate and report different physical properties.",
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
              "results.method.simulation.dft.xc_functional_type": {
                "label": "Jacob's ladder",
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
              "results.method.simulation.precision.apw_cutoff": {
                "label": "APW Cutoff",
                "align": "left"
              },
              "results.method.simulation.precision.basis_set": {
                "label": "Basis Set",
                "align": "left"
              },
              "results.method.simulation.precision.k_line_density": {
                "label": "k-line Density",
                "align": "left"
              },
              "results.method.simulation.precision.native_tier": {
                "label": "Code-specific tier",
                "align": "left"
              },
              "results.method.simulation.precision.planewave_cutoff": {
                "label": "Plane-wave cutoff",
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
              "entry_name": {
                "label": "Name",
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
              "results.method.simulation.dft.xc_functional_type",
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
                "size": "s"
              },
              "elements": {
                "label": "Elements / Formula",
                "level": 1,
                "size": "xl"
              },
              "structure": {
                "label": "Structure / Symmetry",
                "level": 1,
                "size": "s"
              },
              "method": {
                "label": "Method",
                "level": 0,
                "size": "s"
              },
              "precision": {
                "label": "Precision",
                "level": 1,
                "size": "s"
              },
              "dft": {
                "label": "DFT",
                "level": 1,
                "size": "s"
              },
              "gw": {
                "label": "GW",
                "level": 1,
                "size": "s"
              },
              "bse": {
                "label": "BSE",
                "level": 1,
                "size": "s"
              },
              "projection": {
                "label": "Projection",
                "level": 1,
                "size": "s"
              },
              "dmft": {
                "label": "DMFT",
                "level": 1,
                "size": "s"
              },
              "workflow": {
                "label": "Workflow",
                "level": 0,
                "size": "s"
              },
              "molecular_dynamics": {
                "label": "Molecular dynamics",
                "level": 1,
                "size": "s"
              },
              "geometry_optimization": {
                "label": "Geometry Optimization",
                "level": 1,
                "size": "s"
              },
              "properties": {
                "label": "Properties",
                "level": 0,
                "size": "s"
              },
              "electronic": {
                "label": "Electronic",
                "level": 1,
                "size": "s"
              },
              "vibrational": {
                "label": "Vibrational",
                "level": 1,
                "size": "s"
              },
              "mechanical": {
                "label": "Mechanical",
                "level": 1,
                "size": "s"
              },
              "author": {
                "label": "Author / Origin / Dataset",
                "level": 0,
                "size": "m"
              },
              "metadata": {
                "label": "Visibility / IDs / Schema",
                "level": 0,
                "size": "s"
              },
              "optimade": {
                "label": "Optimade",
                "level": 0,
                "size": "m"
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
                    "minH": 3,
                    "minW": 3,
                    "h": 9,
                    "w": 13,
                    "x": 0,
                    "y": 0
                  },
                  "xl": {
                    "minH": 3,
                    "minW": 3,
                    "h": 11,
                    "w": 14,
                    "x": 0,
                    "y": 0
                  },
                  "lg": {
                    "minH": 3,
                    "minW": 3,
                    "h": 11,
                    "w": 14,
                    "x": 0,
                    "y": 0
                  },
                  "md": {
                    "minH": 3,
                    "minW": 3,
                    "h": 8,
                    "w": 12,
                    "x": 0,
                    "y": 0
                  },
                  "sm": {
                    "minH": 3,
                    "minW": 3,
                    "h": 8,
                    "w": 12,
                    "x": 0,
                    "y": 0
                  }
                },
                "quantity": "results.material.elements",
                "scale": "linear"
              },
              {
                "type": "terms",
                "layout": {
                  "xxl": {
                    "minH": 3,
                    "minW": 3,
                    "h": 9,
                    "w": 6,
                    "x": 30,
                    "y": 0
                  },
                  "xl": {
                    "minH": 3,
                    "minW": 3,
                    "h": 6,
                    "w": 6,
                    "x": 24,
                    "y": 5
                  },
                  "lg": {
                    "minH": 3,
                    "minW": 3,
                    "h": 5,
                    "w": 5,
                    "x": 19,
                    "y": 6
                  },
                  "md": {
                    "minH": 3,
                    "minW": 3,
                    "h": 6,
                    "w": 6,
                    "x": 12,
                    "y": 8
                  },
                  "sm": {
                    "minH": 3,
                    "minW": 3,
                    "h": 5,
                    "w": 6,
                    "x": 6,
                    "y": 13
                  }
                },
                "quantity": "results.material.symmetry.space_group_symbol",
                "scale": "linear",
                "showinput": true
              },
              {
                "type": "terms",
                "layout": {
                  "xxl": {
                    "minH": 3,
                    "minW": 3,
                    "h": 9,
                    "w": 6,
                    "x": 19,
                    "y": 0
                  },
                  "xl": {
                    "minH": 3,
                    "minW": 3,
                    "h": 11,
                    "w": 5,
                    "x": 19,
                    "y": 0
                  },
                  "lg": {
                    "minH": 3,
                    "minW": 3,
                    "h": 6,
                    "w": 5,
                    "x": 19,
                    "y": 0
                  },
                  "md": {
                    "minH": 3,
                    "minW": 3,
                    "h": 6,
                    "w": 6,
                    "x": 0,
                    "y": 8
                  },
                  "sm": {
                    "minH": 3,
                    "minW": 3,
                    "h": 5,
                    "w": 6,
                    "x": 6,
                    "y": 8
                  }
                },
                "quantity": "results.material.structural_type",
                "scale": "1/8",
                "showinput": false
              },
              {
                "type": "terms",
                "layout": {
                  "xxl": {
                    "minH": 3,
                    "minW": 3,
                    "h": 9,
                    "w": 6,
                    "x": 13,
                    "y": 0
                  },
                  "xl": {
                    "minH": 3,
                    "minW": 3,
                    "h": 11,
                    "w": 5,
                    "x": 14,
                    "y": 0
                  },
                  "lg": {
                    "minH": 3,
                    "minW": 3,
                    "h": 6,
                    "w": 5,
                    "x": 14,
                    "y": 0
                  },
                  "md": {
                    "minH": 3,
                    "minW": 3,
                    "h": 8,
                    "w": 6,
                    "x": 12,
                    "y": 0
                  },
                  "sm": {
                    "minH": 3,
                    "minW": 3,
                    "h": 5,
                    "w": 6,
                    "x": 0,
                    "y": 8
                  }
                },
                "quantity": "results.method.simulation.program_name",
                "scale": "1/4",
                "showinput": true
              },
              {
                "type": "terms",
                "layout": {
                  "xxl": {
                    "minH": 3,
                    "minW": 3,
                    "h": 9,
                    "w": 5,
                    "x": 25,
                    "y": 0
                  },
                  "xl": {
                    "minH": 3,
                    "minW": 3,
                    "h": 5,
                    "w": 6,
                    "x": 24,
                    "y": 0
                  },
                  "lg": {
                    "minH": 3,
                    "minW": 3,
                    "h": 5,
                    "w": 5,
                    "x": 14,
                    "y": 6
                  },
                  "md": {
                    "minH": 3,
                    "minW": 3,
                    "h": 6,
                    "w": 6,
                    "x": 6,
                    "y": 8
                  },
                  "sm": {
                    "minH": 3,
                    "minW": 3,
                    "h": 5,
                    "w": 6,
                    "x": 0,
                    "y": 13
                  }
                },
                "quantity": "results.material.symmetry.crystal_system",
                "scale": "linear",
                "showinput": false
              }
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
          "category": "Theory",
          "description": "Search materials that are identified from calculations",
          "readme": "This page allows you to search **materials** within\nNOMAD. NOMAD can often automatically detect the material\nfrom individual calculations that contain the full\natomistic structure and can then group the data by using\nthese detected materials. This allows you to search\nindividual materials which have properties that are\naggregated from several entries. Following the link for\na specific material will take you to the corresponding\n[NOMAD Encyclopedia](https://nomad-lab.eu/prod/rae/encyclopedia/#/search)\npage for that material. NOMAD Encyclopedia is a service\nthat is specifically oriented towards materials property\nexploration.\n\nNotice that by default the properties that you search\ncan be combined from several different entries. If\ninstead you wish to search for a material with an\nindividual entry fullfilling your search criteria,\nuncheck the **combine results from several\nentries**-checkbox.",
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
                "size": "s"
              },
              "elements": {
                "label": "Elements / Formula",
                "level": 1,
                "size": "xl"
              },
              "structure": {
                "label": "Structure / Symmetry",
                "level": 1,
                "size": "s"
              },
              "method": {
                "label": "Method",
                "level": 0,
                "size": "s"
              },
              "dft": {
                "label": "DFT",
                "level": 1,
                "size": "s"
              },
              "gw": {
                "label": "GW",
                "level": 1,
                "size": "s"
              },
              "bse": {
                "label": "BSE",
                "level": 1,
                "size": "s"
              },
              "projection": {
                "label": "Projection",
                "level": 1,
                "size": "s"
              },
              "dmft": {
                "label": "DMFT",
                "level": 1,
                "size": "s"
              },
              "workflow": {
                "label": "Workflow",
                "level": 0,
                "size": "s"
              },
              "molecular_dynamics": {
                "label": "Molecular dynamics",
                "level": 1,
                "size": "s"
              },
              "geometry_optimization": {
                "label": "Geometry Optimization",
                "level": 1,
                "size": "s"
              },
              "properties": {
                "label": "Properties",
                "level": 0,
                "size": "s"
              },
              "electronic": {
                "label": "Electronic",
                "level": 1,
                "size": "s"
              },
              "vibrational": {
                "label": "Vibrational",
                "level": 1,
                "size": "s"
              },
              "mechanical": {
                "label": "Mechanical",
                "level": 1,
                "size": "s"
              },
              "author": {
                "label": "Author / Origin / Dataset",
                "level": 0,
                "size": "m"
              },
              "metadata": {
                "label": "Visibility / IDs / Schema",
                "level": 0,
                "size": "s"
              },
              "optimade": {
                "label": "Optimade",
                "level": 0,
                "size": "m"
              },
              "combine": {
                "level": 0,
                "size": "s",
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
          "category": "Experiment",
          "description": "Search electronic lab notebooks",
          "readme": "This page allows you to specifically seach **electronic\nlab notebooks (ELNs)** within NOMAD.  It is very similar\nto the entries search, but with a reduced filter set and\nspecialized arrangement of default columns.",
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
                "size": "s"
              },
              "elements": {
                "label": "Elements / Formula",
                "level": 1,
                "size": "xl"
              },
              "eln": {
                "label": "Electronic Lab Notebook",
                "level": 0,
                "size": "s"
              },
              "custom_quantities": {
                "label": "User Defined Quantities",
                "level": 0,
                "size": "l"
              },
              "author": {
                "label": "Author / Origin / Dataset",
                "level": 0,
                "size": "m"
              },
              "metadata": {
                "label": "Visibility / IDs / Schema",
                "level": 0,
                "size": "s"
              },
              "optimade": {
                "label": "Optimade",
                "level": 0,
                "size": "m"
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
          "category": "Experiment",
          "description": "Search electron energy loss spectroscopy experiments",
          "readme": "This page allows you to spefically search **Electron\nEnergy Loss Spectroscopy (EELS) experiments** within\nNOMAD. It is similar to the entries search, but with a\nreduced filter set and specialized arrangement of\ndefault columns.",
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
              "results.properties.spectroscopic.spectra.provenance.eels.detector_type": {
                "label": "Detector type",
                "align": "left"
              },
              "results.properties.spectroscopic.spectra.provenance.eels.resolution": {
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
              "results.properties.spectroscopic.spectra.provenance.eels.min_energy": {
                "align": "left"
              },
              "results.properties.spectroscopic.spectra.provenance.eels.max_energy": {
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
              "results.properties.spectroscopic.spectra.provenance.eels.detector_type",
              "results.properties.spectroscopic.spectra.provenance.eels.resolution",
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
                "size": "s"
              },
              "elements": {
                "label": "Elements / Formula",
                "level": 1,
                "size": "xl"
              },
              "method": {
                "label": "Method",
                "level": 0,
                "size": "s"
              },
              "eels": {
                "label": "EELS",
                "level": 1,
                "size": "s"
              },
              "author": {
                "label": "Author / Origin / Dataset",
                "level": 0,
                "size": "m"
              },
              "metadata": {
                "label": "Visibility / IDs / Schema",
                "level": 0,
                "size": "s"
              },
              "optimade": {
                "label": "Optimade",
                "level": 0,
                "size": "m"
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
          "category": "Use Cases",
          "description": "Search solar cells",
          "readme": "This page allows you to search **solar cell data**\nwithin NOMAD. The filter menu on the left and the shown\ndefault columns are specifically designed for solar cell\nexploration. The dashboard directly shows useful\ninteractive statistics about the data.",
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
                "size": "s"
              },
              "elements": {
                "label": "Elements / Formula",
                "level": 1,
                "size": "xl"
              },
              "structure": {
                "label": "Structure / Symmetry",
                "level": 1,
                "size": "s"
              },
              "electronic": {
                "label": "Electronic Properties",
                "level": 0,
                "size": "s"
              },
              "solarcell": {
                "label": "Solar Cell Properties",
                "level": 0,
                "size": "s"
              },
              "eln": {
                "label": "Electronic Lab Notebook",
                "level": 0,
                "size": "s"
              },
              "custom_quantities": {
                "label": "User Defined Quantities",
                "level": 0,
                "size": "l"
              },
              "author": {
                "label": "Author / Origin / Dataset",
                "level": 0,
                "size": "m"
              },
              "metadata": {
                "label": "Visibility / IDs / Schema",
                "level": 0,
                "size": "s"
              },
              "optimade": {
                "label": "Optimade",
                "level": 0,
                "size": "m"
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
        },
        "heterogeneouscatalyst": {
          "label": "Heterogeneous Catalysis",
          "path": "heterogeneouscatalyst",
          "resource": "entries",
          "category": "Use Cases",
          "description": "Search heterogeneous catalysts",
          "readme": "This page allows you to search **catalyst and catalysis data**\nwithin NOMAD. The filter menu on the left and the shown\ndefault columns are specifically designed for Heterogeneous Catalyst\nexploration. The dashboard directly shows useful\ninteractive statistics about the data.",
          "pagination": {
            "order_by": "upload_create_time",
            "order": "asc",
            "page_size": 20
          },
          "columns": {
            "options": {
              "results.material.elements": {
                "label": "Elements",
                "align": "left"
              },
              "results.properties.catalytic.catalyst_synthesis.catalyst_type": {
                "label": "Catalyst Type",
                "align": "left"
              },
              "results.properties.catalytic.catalyst_synthesis.preparation_method": {
                "label": "Preparation",
                "align": "left"
              },
              "results.properties.catalytic.catalyst_characterization.surface_area": {
                "label": "Surface Area (m^2/g)",
                "align": "left",
                "format": {
                  "decimals": 2,
                  "mode": "standard"
                }
              },
              "results.properties.catalytic.reactivity.reaction_name": {
                "label": "Reaction Name",
                "align": "left"
              },
              "results.properties.catalytic.reactivity.reaction_class": {
                "label": "Reaction Class",
                "align": "left"
              },
              "results.properties.catalytic.catalyst_synthesis.catalyst_name": {
                "label": "Catalyst Name",
                "align": "left"
              },
              "results.properties.catalytic.reactivity.reactants.name": {
                "label": "Reactants",
                "align": "left"
              },
              "results.properties.catalytic.reactivity.products.name": {
                "label": "Products",
                "align": "left"
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
              "entry_name",
              "results.properties.catalytic.reactivity.reaction_name",
              "results.properties.catalytic.catalyst_synthesis.catalyst_type",
              "results.properties.catalytic.catalyst_synthesis.preparation_method",
              "results.properties.catalytic.catalyst_characterization.surface_area"
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
                "label": "Catalyst Material",
                "level": 0,
                "size": "s"
              },
              "elements": {
                "label": "Elements / Formula",
                "level": 1,
                "size": "xl"
              },
              "structure": {
                "label": "Structure / Symmetry",
                "level": 1,
                "size": "s"
              },
              "heterogeneouscatalyst": {
                "label": "Catalytic Properties",
                "level": 0,
                "size": "s"
              },
              "eln": {
                "label": "Electronic Lab Notebook",
                "level": 0,
                "size": "s"
              },
              "custom_quantities": {
                "label": "User Defined Quantities",
                "level": 0,
                "size": "l"
              },
              "author": {
                "label": "Author / Origin / Dataset",
                "level": 0,
                "size": "m"
              },
              "metadata": {
                "label": "Visibility / IDs / Schema",
                "level": 0,
                "size": "s"
              },
              "optimade": {
                "label": "Optimade",
                "level": 0,
                "size": "m"
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
                    "w": 12,
                    "x": 0,
                    "y": 5
                  },
                  "xl": {
                    "minH": 8,
                    "minW": 12,
                    "h": 8,
                    "w": 12,
                    "x": 0,
                    "y": 5
                  },
                  "lg": {
                    "minH": 8,
                    "minW": 12,
                    "h": 8,
                    "w": 12,
                    "x": 0,
                    "y": 6
                  },
                  "md": {
                    "minH": 8,
                    "minW": 12,
                    "h": 8,
                    "w": 12,
                    "x": 0,
                    "y": 5
                  },
                  "sm": {
                    "minH": 8,
                    "minW": 12,
                    "h": 8,
                    "w": 12,
                    "x": 0,
                    "y": 5
                  }
                },
                "quantity": "results.material.elements",
                "scale": "linear"
              },
              {
                "type": "terms",
                "layout": {
                  "xxl": {
                    "minH": 3,
                    "minW": 3,
                    "h": 5,
                    "w": 6,
                    "x": 6,
                    "y": 0
                  },
                  "xl": {
                    "minH": 3,
                    "minW": 3,
                    "h": 5,
                    "w": 6,
                    "x": 0,
                    "y": 0
                  },
                  "lg": {
                    "minH": 3,
                    "minW": 3,
                    "h": 6,
                    "w": 6,
                    "x": 0,
                    "y": 0
                  },
                  "md": {
                    "minH": 3,
                    "minW": 3,
                    "h": 5,
                    "w": 6,
                    "x": 0,
                    "y": 0
                  },
                  "sm": {
                    "minH": 3,
                    "minW": 3,
                    "h": 5,
                    "w": 4,
                    "x": 0,
                    "y": 0
                  }
                },
                "quantity": "results.properties.catalytic.reactivity.reactants.name",
                "scale": "linear",
                "showinput": true
              },
              {
                "type": "terms",
                "layout": {
                  "xxl": {
                    "minH": 3,
                    "minW": 3,
                    "h": 5,
                    "w": 6,
                    "x": 0,
                    "y": 0
                  },
                  "xl": {
                    "minH": 3,
                    "minW": 3,
                    "h": 5,
                    "w": 6,
                    "x": 12,
                    "y": 0
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
                    "h": 5,
                    "w": 6,
                    "x": 12,
                    "y": 0
                  },
                  "sm": {
                    "minH": 3,
                    "minW": 3,
                    "h": 5,
                    "w": 4,
                    "x": 8,
                    "y": 0
                  }
                },
                "quantity": "results.properties.catalytic.reactivity.reaction_name",
                "scale": "linear",
                "showinput": true
              },
              {
                "type": "terms",
                "layout": {
                  "xxl": {
                    "minH": 3,
                    "minW": 3,
                    "h": 5,
                    "w": 6,
                    "x": 12,
                    "y": 0
                  },
                  "xl": {
                    "minH": 3,
                    "minW": 3,
                    "h": 5,
                    "w": 6,
                    "x": 6,
                    "y": 0
                  },
                  "lg": {
                    "minH": 3,
                    "minW": 3,
                    "h": 6,
                    "w": 6,
                    "x": 6,
                    "y": 0
                  },
                  "md": {
                    "minH": 3,
                    "minW": 3,
                    "h": 5,
                    "w": 6,
                    "x": 6,
                    "y": 0
                  },
                  "sm": {
                    "minH": 3,
                    "minW": 3,
                    "h": 5,
                    "w": 4,
                    "x": 4,
                    "y": 0
                  }
                },
                "quantity": "results.properties.catalytic.reactivity.products.name",
                "scale": "linear",
                "showinput": true
              },
              {
                "type": "terms",
                "layout": {
                  "xxl": {
                    "minH": 3,
                    "minW": 3,
                    "h": 4,
                    "w": 6,
                    "x": 12,
                    "y": 5
                  },
                  "xl": {
                    "minH": 3,
                    "minW": 3,
                    "h": 4,
                    "w": 6,
                    "x": 12,
                    "y": 5
                  },
                  "lg": {
                    "minH": 3,
                    "minW": 3,
                    "h": 4,
                    "w": 6,
                    "x": 12,
                    "y": 6
                  },
                  "md": {
                    "minH": 3,
                    "minW": 3,
                    "h": 4,
                    "w": 6,
                    "x": 12,
                    "y": 5
                  },
                  "sm": {
                    "minH": 3,
                    "minW": 3,
                    "h": 3,
                    "w": 4,
                    "x": 8,
                    "y": 13
                  }
                },
                "quantity": "results.properties.catalytic.catalyst_synthesis.preparation_method",
                "scale": "linear",
                "showinput": true
              },
              {
                "type": "terms",
                "layout": {
                  "xxl": {
                    "minH": 3,
                    "minW": 3,
                    "h": 4,
                    "w": 6,
                    "x": 12,
                    "y": 9
                  },
                  "xl": {
                    "minH": 3,
                    "minW": 3,
                    "h": 4,
                    "w": 6,
                    "x": 12,
                    "y": 9
                  },
                  "lg": {
                    "minH": 3,
                    "minW": 3,
                    "h": 4,
                    "w": 6,
                    "x": 12,
                    "y": 10
                  },
                  "md": {
                    "minH": 3,
                    "minW": 3,
                    "h": 4,
                    "w": 6,
                    "x": 12,
                    "y": 9
                  },
                  "sm": {
                    "minH": 3,
                    "minW": 3,
                    "h": 3,
                    "w": 4,
                    "x": 8,
                    "y": 16
                  }
                },
                "quantity": "results.properties.catalytic.catalyst_synthesis.catalyst_type",
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
                    "w": 9,
                    "x": 0,
                    "y": 13
                  },
                  "xl": {
                    "minH": 3,
                    "minW": 8,
                    "h": 4,
                    "w": 9,
                    "x": 0,
                    "y": 13
                  },
                  "lg": {
                    "minH": 3,
                    "minW": 8,
                    "h": 4,
                    "w": 9,
                    "x": 0,
                    "y": 14
                  },
                  "md": {
                    "minH": 3,
                    "minW": 8,
                    "h": 3,
                    "w": 9,
                    "x": 0,
                    "y": 13
                  },
                  "sm": {
                    "minH": 3,
                    "minW": 8,
                    "h": 3,
                    "w": 8,
                    "x": 0,
                    "y": 13
                  }
                },
                "quantity": "results.properties.catalytic.reactivity.test_temperatures",
                "scale": "linear",
                "autorange": false,
                "showinput": false,
                "nbins": 30
              },
              {
                "type": "histogram",
                "layout": {
                  "xxl": {
                    "minH": 3,
                    "minW": 8,
                    "h": 3,
                    "w": 9,
                    "x": 0,
                    "y": 16
                  },
                  "xl": {
                    "minH": 3,
                    "minW": 8,
                    "h": 4,
                    "w": 9,
                    "x": 0,
                    "y": 17
                  },
                  "lg": {
                    "minH": 3,
                    "minW": 8,
                    "h": 4,
                    "w": 9,
                    "x": 0,
                    "y": 18
                  },
                  "md": {
                    "minH": 3,
                    "minW": 8,
                    "h": 3,
                    "w": 9,
                    "x": 9,
                    "y": 16
                  },
                  "sm": {
                    "minH": 3,
                    "minW": 8,
                    "h": 3,
                    "w": 8,
                    "x": 0,
                    "y": 22
                  }
                },
                "quantity": "results.properties.catalytic.reactivity.gas_hourly_space_velocity",
                "scale": "linear",
                "autorange": false,
                "showinput": false,
                "nbins": 30
              },
              {
                "type": "histogram",
                "layout": {
                  "xxl": {
                    "minH": 3,
                    "minW": 8,
                    "h": 3,
                    "w": 9,
                    "x": 9,
                    "y": 13
                  },
                  "xl": {
                    "minH": 3,
                    "minW": 8,
                    "h": 4,
                    "w": 9,
                    "x": 9,
                    "y": 13
                  },
                  "lg": {
                    "minH": 3,
                    "minW": 8,
                    "h": 4,
                    "w": 9,
                    "x": 9,
                    "y": 14
                  },
                  "md": {
                    "minH": 3,
                    "minW": 8,
                    "h": 3,
                    "w": 9,
                    "x": 9,
                    "y": 13
                  },
                  "sm": {
                    "minH": 3,
                    "minW": 8,
                    "h": 3,
                    "w": 8,
                    "x": 0,
                    "y": 16
                  }
                },
                "quantity": "results.properties.catalytic.reactivity.reactants.gas_concentration_in",
                "scale": "linear",
                "autorange": false,
                "showinput": false,
                "nbins": 30
              },
              {
                "type": "histogram",
                "layout": {
                  "xxl": {
                    "minH": 3,
                    "minW": 8,
                    "h": 3,
                    "w": 9,
                    "x": 9,
                    "y": 16
                  },
                  "xl": {
                    "minH": 3,
                    "minW": 8,
                    "h": 4,
                    "w": 9,
                    "x": 9,
                    "y": 17
                  },
                  "lg": {
                    "minH": 3,
                    "minW": 8,
                    "h": 4,
                    "w": 9,
                    "x": 9,
                    "y": 14
                  },
                  "md": {
                    "minH": 3,
                    "minW": 8,
                    "h": 3,
                    "w": 9,
                    "x": 0,
                    "y": 16
                  },
                  "sm": {
                    "minH": 3,
                    "minW": 8,
                    "h": 3,
                    "w": 8,
                    "x": 0,
                    "y": 16
                  }
                },
                "quantity": "results.properties.catalytic.reactivity.pressure",
                "scale": "linear",
                "autorange": false,
                "showinput": false,
                "nbins": 30
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
                    "y": 19
                  },
                  "xl": {
                    "minH": 3,
                    "minW": 8,
                    "h": 4,
                    "w": 9,
                    "x": 0,
                    "y": 21
                  },
                  "lg": {
                    "minH": 3,
                    "minW": 8,
                    "h": 4,
                    "w": 9,
                    "x": 0,
                    "y": 26
                  },
                  "md": {
                    "minH": 3,
                    "minW": 8,
                    "h": 3,
                    "w": 9,
                    "x": 0,
                    "y": 22
                  },
                  "sm": {
                    "minH": 3,
                    "minW": 8,
                    "h": 3,
                    "w": 8,
                    "x": 0,
                    "y": 33
                  }
                },
                "quantity": "results.properties.catalytic.reactivity.products.selectivity",
                "scale": "linear",
                "autorange": false,
                "showinput": false,
                "nbins": 30
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
                    "y": 22
                  },
                  "xl": {
                    "minH": 3,
                    "minW": 8,
                    "h": 4,
                    "w": 9,
                    "x": 0,
                    "y": 25
                  },
                  "lg": {
                    "minH": 3,
                    "minW": 8,
                    "h": 4,
                    "w": 9,
                    "x": 0,
                    "y": 22
                  },
                  "md": {
                    "minH": 3,
                    "minW": 8,
                    "h": 3,
                    "w": 9,
                    "x": 0,
                    "y": 19
                  },
                  "sm": {
                    "minH": 3,
                    "minW": 8,
                    "h": 3,
                    "w": 8,
                    "x": 0,
                    "y": 30
                  }
                },
                "quantity": "results.properties.catalytic.reactivity.reactants.conversion",
                "scale": "linear",
                "autorange": false,
                "showinput": false,
                "nbins": 30
              },
              {
                "type": "histogram",
                "layout": {
                  "xxl": {
                    "minH": 3,
                    "minW": 8,
                    "h": 3,
                    "w": 8,
                    "x": 8,
                    "y": 25
                  },
                  "xl": {
                    "minH": 3,
                    "minW": 8,
                    "h": 4,
                    "w": 9,
                    "x": 9,
                    "y": 29
                  },
                  "lg": {
                    "minH": 3,
                    "minW": 8,
                    "h": 4,
                    "w": 12,
                    "x": 0,
                    "y": 30
                  },
                  "md": {
                    "minH": 3,
                    "minW": 8,
                    "h": 3,
                    "w": 9,
                    "x": 0,
                    "y": 25
                  },
                  "sm": {
                    "minH": 3,
                    "minW": 8,
                    "h": 3,
                    "w": 8,
                    "x": 0,
                    "y": 36
                  }
                },
                "quantity": "results.properties.catalytic.reactivity.rates.reaction_rate",
                "scale": "linear",
                "autorange": false,
                "showinput": false,
                "nbins": 30
              },
              {
                "type": "scatterplot",
                "layout": {
                  "xxl": {
                    "minH": 3,
                    "minW": 3,
                    "h": 6,
                    "w": 10,
                    "x": 8,
                    "y": 19
                  },
                  "xl": {
                    "minH": 3,
                    "minW": 3,
                    "h": 8,
                    "w": 9,
                    "x": 9,
                    "y": 21
                  },
                  "lg": {
                    "minH": 3,
                    "minW": 3,
                    "h": 8,
                    "w": 9,
                    "x": 9,
                    "y": 22
                  },
                  "md": {
                    "minH": 3,
                    "minW": 3,
                    "h": 6,
                    "w": 9,
                    "x": 9,
                    "y": 19
                  },
                  "sm": {
                    "minH": 3,
                    "minW": 3,
                    "h": 5,
                    "w": 8,
                    "x": 9,
                    "y": 25
                  }
                },
                "x": "results.properties.catalytic.reactivity.reactants.conversion",
                "y": "results.properties.catalytic.reactivity.products.selectivity",
                "color": "results.properties.catalytic.catalyst_characterization.surface_area",
                "size": 1000,
                "autorange": true
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
                    "y": 25
                  },
                  "xl": {
                    "minH": 3,
                    "minW": 8,
                    "h": 4,
                    "w": 9,
                    "x": 0,
                    "y": 29
                  },
                  "lg": {
                    "minH": 3,
                    "minW": 8,
                    "h": 4,
                    "w": 12,
                    "x": 0,
                    "y": 34
                  },
                  "md": {
                    "minH": 3,
                    "minW": 8,
                    "h": 3,
                    "w": 9,
                    "x": 0,
                    "y": 28
                  },
                  "sm": {
                    "minH": 3,
                    "minW": 8,
                    "h": 3,
                    "w": 8,
                    "x": 0,
                    "y": 39
                  }
                },
                "quantity": "results.properties.catalytic.catalyst_characterization.surface_area",
                "scale": "1/4",
                "autorange": false,
                "showinput": false,
                "nbins": 30
              }
            ]
          },
          "filters_locked": {
            "quantities": "results.properties.catalytic"
          }
        }
      }
    },
    "north": {},
    "example_uploads": {}
  }
}
