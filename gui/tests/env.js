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
      "options": {
        "Custom": {
          "label": "Custom",
          "units": {
            "angle": {
              "definition": "\u00b0",
              "locked": false
            },
            "energy": {
              "definition": "eV",
              "locked": false
            },
            "length": {
              "definition": "\u00c5",
              "locked": false
            },
            "pressure": {
              "definition": "GPa",
              "locked": false
            },
            "time": {
              "definition": "fs",
              "locked": false
            },
            "dimensionless": {
              "definition": "dimensionless",
              "locked": false
            },
            "mass": {
              "definition": "kg",
              "locked": false
            },
            "current": {
              "definition": "A",
              "locked": false
            },
            "temperature": {
              "definition": "K",
              "locked": false
            },
            "luminosity": {
              "definition": "cd",
              "locked": false
            },
            "luminous_flux": {
              "definition": "lm",
              "locked": false
            },
            "substance": {
              "definition": "mol",
              "locked": false
            },
            "information": {
              "definition": "bit",
              "locked": false
            },
            "force": {
              "definition": "N",
              "locked": false
            },
            "power": {
              "definition": "W",
              "locked": false
            },
            "charge": {
              "definition": "C",
              "locked": false
            },
            "resistance": {
              "definition": "\u03a9",
              "locked": false
            },
            "conductance": {
              "definition": "S",
              "locked": false
            },
            "inductance": {
              "definition": "H",
              "locked": false
            },
            "magnetic_flux": {
              "definition": "Wb",
              "locked": false
            },
            "magnetic_field": {
              "definition": "T",
              "locked": false
            },
            "frequency": {
              "definition": "Hz",
              "locked": false
            },
            "luminance": {
              "definition": "nit",
              "locked": false
            },
            "illuminance": {
              "definition": "lx",
              "locked": false
            },
            "electric_potential": {
              "definition": "V",
              "locked": false
            },
            "capacitance": {
              "definition": "F",
              "locked": false
            },
            "activity": {
              "definition": "kat",
              "locked": false
            }
          }
        },
        "SI": {
          "label": "International System of Units (SI)",
          "units": {
            "activity": {
              "definition": "kat",
              "locked": true
            },
            "angle": {
              "definition": "rad",
              "locked": true
            },
            "capacitance": {
              "definition": "F",
              "locked": true
            },
            "charge": {
              "definition": "C",
              "locked": true
            },
            "conductance": {
              "definition": "S",
              "locked": true
            },
            "current": {
              "definition": "A",
              "locked": true
            },
            "dimensionless": {
              "definition": "dimensionless",
              "locked": true
            },
            "electric_potential": {
              "definition": "V",
              "locked": true
            },
            "energy": {
              "definition": "J",
              "locked": true
            },
            "force": {
              "definition": "N",
              "locked": true
            },
            "frequency": {
              "definition": "Hz",
              "locked": true
            },
            "illuminance": {
              "definition": "lx",
              "locked": true
            },
            "inductance": {
              "definition": "H",
              "locked": true
            },
            "information": {
              "definition": "bit",
              "locked": true
            },
            "length": {
              "definition": "m",
              "locked": true
            },
            "luminance": {
              "definition": "nit",
              "locked": true
            },
            "luminosity": {
              "definition": "cd",
              "locked": true
            },
            "luminous_flux": {
              "definition": "lm",
              "locked": true
            },
            "magnetic_field": {
              "definition": "T",
              "locked": true
            },
            "magnetic_flux": {
              "definition": "Wb",
              "locked": true
            },
            "mass": {
              "definition": "kg",
              "locked": true
            },
            "power": {
              "definition": "W",
              "locked": true
            },
            "pressure": {
              "definition": "Pa",
              "locked": true
            },
            "resistance": {
              "definition": "\u03a9",
              "locked": true
            },
            "substance": {
              "definition": "mol",
              "locked": true
            },
            "temperature": {
              "definition": "K",
              "locked": true
            },
            "time": {
              "definition": "s",
              "locked": true
            }
          }
        },
        "AU": {
          "label": "Hartree atomic units (AU)",
          "units": {
            "activity": {
              "definition": "kat",
              "locked": false
            },
            "angle": {
              "definition": "rad",
              "locked": false
            },
            "capacitance": {
              "definition": "F",
              "locked": false
            },
            "charge": {
              "definition": "C",
              "locked": false
            },
            "conductance": {
              "definition": "S",
              "locked": false
            },
            "current": {
              "definition": "atomic_unit_of_current",
              "locked": true
            },
            "dimensionless": {
              "definition": "dimensionless",
              "locked": true
            },
            "electric_potential": {
              "definition": "V",
              "locked": false
            },
            "energy": {
              "definition": "Ha",
              "locked": true
            },
            "force": {
              "definition": "atomic_unit_of_force",
              "locked": true
            },
            "frequency": {
              "definition": "Hz",
              "locked": false
            },
            "illuminance": {
              "definition": "lx",
              "locked": false
            },
            "inductance": {
              "definition": "H",
              "locked": false
            },
            "information": {
              "definition": "bit",
              "locked": false
            },
            "length": {
              "definition": "bohr",
              "locked": true
            },
            "luminance": {
              "definition": "nit",
              "locked": false
            },
            "luminosity": {
              "definition": "cd",
              "locked": false
            },
            "luminous_flux": {
              "definition": "lm",
              "locked": false
            },
            "magnetic_field": {
              "definition": "T",
              "locked": false
            },
            "magnetic_flux": {
              "definition": "Wb",
              "locked": false
            },
            "mass": {
              "definition": "m_e",
              "locked": true
            },
            "power": {
              "definition": "W",
              "locked": false
            },
            "pressure": {
              "definition": "atomic_unit_of_pressure",
              "locked": true
            },
            "resistance": {
              "definition": "\u03a9",
              "locked": false
            },
            "substance": {
              "definition": "mol",
              "locked": false
            },
            "temperature": {
              "definition": "atomic_unit_of_temperature",
              "locked": true
            },
            "time": {
              "definition": "atomic_unit_of_time",
              "locked": true
            }
          }
        }
      },
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
          "history": {
            "error": "Could not render history card."
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
          "readme": "This page allows you to search **entries** within NOMAD. Entries represent any individual data items that have been uploaded to NOMAD, no matter whether they come from theoretical calculations, experiments, lab notebooks or any other source of data. This allows you to perform cross-domain queries, but if you are interested in a specific subfield, you should see if a specific application exists for it in the explore menu to get more details.",
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
                "align": "left"
              },
              "entry_create_time": {
                "align": "left"
              },
              "upload_name": {
                "align": "left"
              },
              "upload_id": {
                "align": "left"
              },
              "upload_create_time": {
                "label": "Upload time",
                "align": "left"
              },
              "authors": {
                "align": "left"
              },
              "results.method.method_name": {
                "align": "left"
              },
              "results.method.simulation.program_name": {
                "align": "left"
              },
              "results.method.simulation.dft.xc_functional_type": {
                "align": "left"
              },
              "results.method.simulation.precision.apw_cutoff": {
                "align": "left"
              },
              "results.method.simulation.precision.basis_set": {
                "align": "left"
              },
              "results.method.simulation.precision.k_line_density": {
                "align": "left"
              },
              "results.method.simulation.precision.native_tier": {
                "align": "left"
              },
              "results.method.simulation.precision.planewave_cutoff": {
                "align": "left"
              },
              "results.material.structural_type": {
                "align": "left"
              },
              "results.material.symmetry.crystal_system": {
                "align": "left"
              },
              "results.material.symmetry.space_group_symbol": {
                "align": "left"
              },
              "results.material.symmetry.space_group_number": {
                "align": "left"
              },
              "results.eln.lab_ids": {
                "align": "left"
              },
              "results.eln.sections": {
                "align": "left"
              },
              "results.eln.methods": {
                "align": "left"
              },
              "results.eln.tags": {
                "align": "left"
              },
              "results.eln.instruments": {
                "align": "left"
              },
              "mainfile": {
                "align": "left"
              },
              "comment": {
                "align": "left"
              },
              "references": {
                "align": "left"
              },
              "datasets": {
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
              "tb": {
                "label": "TB",
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
          },
          "search_syntaxes": {
            "exclude": [
              "free_text"
            ]
          }
        },
        "calculations": {
          "label": "Calculations",
          "path": "calculations",
          "resource": "entries",
          "category": "Theory",
          "description": "Search calculations",
          "readme": "This page allows you to search **calculations** within NOMAD. Calculations typically come from a specific simulation software that uses an approximate model to investigate and report different physical properties.",
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
                "align": "left"
              },
              "results.method.method_name": {
                "align": "left"
              },
              "results.method.simulation.dft.xc_functional_type": {
                "align": "left"
              },
              "upload_create_time": {
                "label": "Upload time",
                "align": "left"
              },
              "authors": {
                "align": "left"
              },
              "results.method.simulation.precision.apw_cutoff": {
                "align": "left"
              },
              "results.method.simulation.precision.basis_set": {
                "align": "left"
              },
              "results.method.simulation.precision.k_line_density": {
                "align": "left"
              },
              "results.method.simulation.precision.native_tier": {
                "align": "left"
              },
              "results.method.simulation.precision.planewave_cutoff": {
                "align": "left"
              },
              "results.material.structural_type": {
                "align": "left"
              },
              "results.material.symmetry.crystal_system": {
                "align": "left"
              },
              "results.material.symmetry.space_group_symbol": {
                "align": "left"
              },
              "results.material.symmetry.space_group_number": {
                "align": "left"
              },
              "entry_name": {
                "label": "Name",
                "align": "left"
              },
              "mainfile": {
                "align": "left"
              },
              "comment": {
                "align": "left"
              },
              "references": {
                "align": "left"
              },
              "datasets": {
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
              "tb": {
                "label": "TB",
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
                  "lg": {
                    "h": 11,
                    "w": 14,
                    "x": 0,
                    "y": 0,
                    "minH": 3,
                    "minW": 3
                  },
                  "md": {
                    "h": 8,
                    "w": 12,
                    "x": 0,
                    "y": 0,
                    "minH": 3,
                    "minW": 3
                  },
                  "sm": {
                    "h": 8,
                    "w": 12,
                    "x": 0,
                    "y": 0,
                    "minH": 3,
                    "minW": 3
                  },
                  "xl": {
                    "h": 11,
                    "w": 14,
                    "x": 0,
                    "y": 0,
                    "minH": 3,
                    "minW": 3
                  },
                  "xxl": {
                    "h": 9,
                    "w": 13,
                    "x": 0,
                    "y": 0,
                    "minH": 3,
                    "minW": 3
                  }
                },
                "quantity": "results.material.elements",
                "scale": "linear"
              },
              {
                "type": "terms",
                "layout": {
                  "lg": {
                    "h": 5,
                    "w": 5,
                    "x": 19,
                    "y": 6,
                    "minH": 3,
                    "minW": 3
                  },
                  "md": {
                    "h": 6,
                    "w": 6,
                    "x": 12,
                    "y": 8,
                    "minH": 3,
                    "minW": 3
                  },
                  "sm": {
                    "h": 5,
                    "w": 6,
                    "x": 6,
                    "y": 13,
                    "minH": 3,
                    "minW": 3
                  },
                  "xl": {
                    "h": 6,
                    "w": 6,
                    "x": 24,
                    "y": 5,
                    "minH": 3,
                    "minW": 3
                  },
                  "xxl": {
                    "h": 9,
                    "w": 6,
                    "x": 30,
                    "y": 0,
                    "minH": 3,
                    "minW": 3
                  }
                },
                "quantity": "results.material.symmetry.space_group_symbol",
                "scale": "linear",
                "showinput": true
              },
              {
                "type": "terms",
                "layout": {
                  "lg": {
                    "h": 6,
                    "w": 5,
                    "x": 19,
                    "y": 0,
                    "minH": 3,
                    "minW": 3
                  },
                  "md": {
                    "h": 6,
                    "w": 6,
                    "x": 0,
                    "y": 8,
                    "minH": 3,
                    "minW": 3
                  },
                  "sm": {
                    "h": 5,
                    "w": 6,
                    "x": 6,
                    "y": 8,
                    "minH": 3,
                    "minW": 3
                  },
                  "xl": {
                    "h": 11,
                    "w": 5,
                    "x": 19,
                    "y": 0,
                    "minH": 3,
                    "minW": 3
                  },
                  "xxl": {
                    "h": 9,
                    "w": 6,
                    "x": 19,
                    "y": 0,
                    "minH": 3,
                    "minW": 3
                  }
                },
                "quantity": "results.material.structural_type",
                "scale": "1/8",
                "showinput": false
              },
              {
                "type": "terms",
                "layout": {
                  "lg": {
                    "h": 6,
                    "w": 5,
                    "x": 14,
                    "y": 0,
                    "minH": 3,
                    "minW": 3
                  },
                  "md": {
                    "h": 8,
                    "w": 6,
                    "x": 12,
                    "y": 0,
                    "minH": 3,
                    "minW": 3
                  },
                  "sm": {
                    "h": 5,
                    "w": 6,
                    "x": 0,
                    "y": 8,
                    "minH": 3,
                    "minW": 3
                  },
                  "xl": {
                    "h": 11,
                    "w": 5,
                    "x": 14,
                    "y": 0,
                    "minH": 3,
                    "minW": 3
                  },
                  "xxl": {
                    "h": 9,
                    "w": 6,
                    "x": 13,
                    "y": 0,
                    "minH": 3,
                    "minW": 3
                  }
                },
                "quantity": "results.method.simulation.program_name",
                "scale": "1/4",
                "showinput": true
              },
              {
                "type": "terms",
                "layout": {
                  "lg": {
                    "h": 5,
                    "w": 5,
                    "x": 14,
                    "y": 6,
                    "minH": 3,
                    "minW": 3
                  },
                  "md": {
                    "h": 6,
                    "w": 6,
                    "x": 6,
                    "y": 8,
                    "minH": 3,
                    "minW": 3
                  },
                  "sm": {
                    "h": 5,
                    "w": 6,
                    "x": 0,
                    "y": 13,
                    "minH": 3,
                    "minW": 3
                  },
                  "xl": {
                    "h": 5,
                    "w": 6,
                    "x": 24,
                    "y": 0,
                    "minH": 3,
                    "minW": 3
                  },
                  "xxl": {
                    "h": 9,
                    "w": 5,
                    "x": 25,
                    "y": 0,
                    "minH": 3,
                    "minW": 3
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
          },
          "search_syntaxes": {
            "exclude": [
              "free_text"
            ]
          }
        },
        "materials": {
          "label": "Materials",
          "path": "materials",
          "resource": "materials",
          "category": "Theory",
          "description": "Search materials that are identified from calculations",
          "readme": "This page allows you to search **materials** within NOMAD. NOMAD can often automatically detect the material from individual calculations that contain the full atomistic structure and can then group the data by using these detected materials. This allows you to search individual materials which have properties that are aggregated from several entries. Following the link for a specific material will take you to the corresponding [NOMAD Encyclopedia](https://nomad-lab.eu/prod/rae/encyclopedia/#/search) page for that material. NOMAD Encyclopedia is a service that is specifically oriented towards materials property exploration.\nNotice that by default the properties that you search can be combined from several different entries. If instead you wish to search for a material with an individual entry fullfilling your search criteria, uncheck the **combine results from several entries**-checkbox.",
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
                "align": "left"
              },
              "symmetry.structure_name": {
                "align": "left"
              },
              "symmetry.space_group_number": {
                "align": "left"
              },
              "symmetry.crystal_system": {
                "align": "left"
              },
              "symmetry.space_group_symbol": {
                "align": "left"
              },
              "material_id": {
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
              "tb": {
                "label": "TB",
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
          },
          "search_syntaxes": {
            "exclude": [
              "free_text"
            ]
          }
        },
        "eln": {
          "label": "ELN",
          "path": "eln",
          "resource": "entries",
          "category": "Experiment",
          "description": "Search electronic lab notebooks",
          "readme": "This page allows you to specifically seach **electronic lab notebooks (ELNs)** within NOMAD. It is very similar to the entries search, but with a reduced filter set and specialized arrangement of default columns.",
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
                "align": "left"
              },
              "upload_create_time": {
                "label": "Upload time",
                "align": "left"
              },
              "authors": {
                "align": "left"
              },
              "results.material.chemical_formula_hill": {
                "label": "Formula",
                "align": "left"
              },
              "results.method.method_name": {
                "align": "left"
              },
              "results.eln.lab_ids": {
                "align": "left"
              },
              "results.eln.sections": {
                "align": "left"
              },
              "results.eln.methods": {
                "align": "left"
              },
              "results.eln.tags": {
                "align": "left"
              },
              "results.eln.instruments": {
                "align": "left"
              },
              "mainfile": {
                "align": "left"
              },
              "comment": {
                "align": "left"
              },
              "references": {
                "align": "left"
              },
              "datasets": {
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
          "readme": "This page allows you to spefically search **Electron Energy Loss Spectroscopy (EELS) experiments** within NOMAD. It is similar to the entries search, but with a reduced filter set and specialized arrangement of default columns.",
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
                "align": "left"
              },
              "results.properties.spectroscopic.spectra.provenance.eels.resolution": {
                "align": "left"
              },
              "upload_create_time": {
                "label": "Upload time",
                "align": "left"
              },
              "authors": {
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
                "align": "left"
              },
              "mainfile": {
                "align": "left"
              },
              "comment": {
                "align": "left"
              },
              "references": {
                "align": "left"
              },
              "datasets": {
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
          },
          "search_syntaxes": {
            "exclude": [
              "free_text"
            ]
          }
        },
        "solarcells": {
          "label": "Solar Cells",
          "path": "solarcells",
          "resource": "entries",
          "category": "Use Cases",
          "description": "Search solar cells",
          "readme": "This page allows you to search **solar cell data** within NOMAD. The filter menu on the left and the shown default columns are specifically designed for solar cell exploration. The dashboard directly shows useful interactive statistics about the data.",
          "pagination": {
            "order_by": "results.properties.optoelectronic.solar_cell.efficiency",
            "order": "desc",
            "page_size": 20
          },
          "columns": {
            "options": {
              "results.material.chemical_formula_descriptive": {
                "label": "Descriptive formula",
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
                "align": "left",
                "unit": "V",
                "format": {
                  "decimals": 3,
                  "mode": "standard"
                }
              },
              "results.properties.optoelectronic.solar_cell.short_circuit_current_density": {
                "align": "left",
                "unit": "A/m**2",
                "format": {
                  "decimals": 3,
                  "mode": "standard"
                }
              },
              "results.properties.optoelectronic.solar_cell.fill_factor": {
                "align": "left",
                "format": {
                  "decimals": 3,
                  "mode": "standard"
                }
              },
              "references": {
                "align": "left"
              },
              "results.material.chemical_formula_hill": {
                "label": "Formula",
                "align": "left"
              },
              "results.material.structural_type": {
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
                "align": "left"
              },
              "results.eln.sections": {
                "align": "left"
              },
              "results.eln.methods": {
                "align": "left"
              },
              "results.eln.tags": {
                "align": "left"
              },
              "results.eln.instruments": {
                "align": "left"
              },
              "entry_name": {
                "label": "Name",
                "align": "left"
              },
              "entry_type": {
                "align": "left"
              },
              "mainfile": {
                "align": "left"
              },
              "upload_create_time": {
                "label": "Upload time",
                "align": "left"
              },
              "authors": {
                "align": "left"
              },
              "comment": {
                "align": "left"
              },
              "datasets": {
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
            "include": [
              "*#perovskite_solar_cell_database.schema.PerovskiteSolarCell"
            ],
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
                  "lg": {
                    "h": 8,
                    "w": 12,
                    "x": 0,
                    "y": 0,
                    "minH": 8,
                    "minW": 12
                  },
                  "md": {
                    "h": 8,
                    "w": 12,
                    "x": 0,
                    "y": 0,
                    "minH": 8,
                    "minW": 12
                  },
                  "sm": {
                    "h": 8,
                    "w": 12,
                    "x": 0,
                    "y": 16,
                    "minH": 8,
                    "minW": 12
                  },
                  "xl": {
                    "h": 8,
                    "w": 12,
                    "x": 0,
                    "y": 0,
                    "minH": 8,
                    "minW": 12
                  },
                  "xxl": {
                    "h": 8,
                    "w": 13,
                    "x": 0,
                    "y": 0,
                    "minH": 8,
                    "minW": 12
                  }
                },
                "quantity": "results.material.elements",
                "scale": "linear"
              },
              {
                "type": "scatterplot",
                "layout": {
                  "lg": {
                    "h": 6,
                    "w": 12,
                    "x": 0,
                    "y": 8,
                    "minH": 3,
                    "minW": 3
                  },
                  "md": {
                    "h": 6,
                    "w": 9,
                    "x": 0,
                    "y": 8,
                    "minH": 3,
                    "minW": 3
                  },
                  "sm": {
                    "h": 5,
                    "w": 6,
                    "x": 0,
                    "y": 0,
                    "minH": 3,
                    "minW": 3
                  },
                  "xl": {
                    "h": 8,
                    "w": 9,
                    "x": 12,
                    "y": 0,
                    "minH": 3,
                    "minW": 3
                  },
                  "xxl": {
                    "h": 8,
                    "w": 12,
                    "x": 24,
                    "y": 0,
                    "minH": 3,
                    "minW": 3
                  }
                },
                "x": {
                  "quantity": "results.properties.optoelectronic.solar_cell.open_circuit_voltage"
                },
                "y": {
                  "title": "Efficiency (%)",
                  "quantity": "results.properties.optoelectronic.solar_cell.efficiency"
                },
                "markers": {
                  "color": {
                    "unit": "mA/cm^2",
                    "quantity": "results.properties.optoelectronic.solar_cell.short_circuit_current_density"
                  }
                },
                "size": 1000,
                "autorange": true
              },
              {
                "type": "scatterplot",
                "layout": {
                  "lg": {
                    "h": 6,
                    "w": 12,
                    "x": 0,
                    "y": 14,
                    "minH": 3,
                    "minW": 3
                  },
                  "md": {
                    "h": 6,
                    "w": 9,
                    "x": 9,
                    "y": 8,
                    "minH": 3,
                    "minW": 3
                  },
                  "sm": {
                    "h": 5,
                    "w": 6,
                    "x": 6,
                    "y": 0,
                    "minH": 3,
                    "minW": 3
                  },
                  "xl": {
                    "h": 8,
                    "w": 9,
                    "x": 21,
                    "y": 0,
                    "minH": 3,
                    "minW": 3
                  },
                  "xxl": {
                    "h": 8,
                    "w": 11,
                    "x": 13,
                    "y": 0,
                    "minH": 3,
                    "minW": 3
                  }
                },
                "x": {
                  "quantity": "results.properties.optoelectronic.solar_cell.open_circuit_voltage"
                },
                "y": {
                  "title": "Efficiency (%)",
                  "quantity": "results.properties.optoelectronic.solar_cell.efficiency"
                },
                "markers": {
                  "color": {
                    "quantity": "results.properties.optoelectronic.solar_cell.device_architecture"
                  }
                },
                "size": 1000,
                "autorange": true
              },
              {
                "type": "terms",
                "layout": {
                  "lg": {
                    "h": 6,
                    "w": 6,
                    "x": 12,
                    "y": 0,
                    "minH": 3,
                    "minW": 3
                  },
                  "md": {
                    "h": 4,
                    "w": 6,
                    "x": 12,
                    "y": 4,
                    "minH": 3,
                    "minW": 3
                  },
                  "sm": {
                    "h": 6,
                    "w": 4,
                    "x": 0,
                    "y": 10,
                    "minH": 3,
                    "minW": 3
                  },
                  "xl": {
                    "h": 6,
                    "w": 6,
                    "x": 14,
                    "y": 8,
                    "minH": 3,
                    "minW": 3
                  },
                  "xxl": {
                    "h": 6,
                    "w": 6,
                    "x": 14,
                    "y": 8,
                    "minH": 3,
                    "minW": 3
                  }
                },
                "quantity": "results.properties.optoelectronic.solar_cell.device_stack",
                "scale": "linear",
                "showinput": true
              },
              {
                "type": "histogram",
                "layout": {
                  "lg": {
                    "h": 4,
                    "w": 12,
                    "x": 12,
                    "y": 12,
                    "minH": 3,
                    "minW": 3
                  },
                  "md": {
                    "h": 3,
                    "w": 8,
                    "x": 10,
                    "y": 17,
                    "minH": 3,
                    "minW": 3
                  },
                  "sm": {
                    "h": 3,
                    "w": 8,
                    "x": 4,
                    "y": 13,
                    "minH": 3,
                    "minW": 3
                  },
                  "xl": {
                    "h": 3,
                    "w": 8,
                    "x": 0,
                    "y": 11,
                    "minH": 3,
                    "minW": 3
                  },
                  "xxl": {
                    "h": 3,
                    "w": 8,
                    "x": 0,
                    "y": 8,
                    "minH": 3,
                    "minW": 3
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
                  "lg": {
                    "h": 6,
                    "w": 6,
                    "x": 18,
                    "y": 0,
                    "minH": 3,
                    "minW": 3
                  },
                  "md": {
                    "h": 4,
                    "w": 6,
                    "x": 12,
                    "y": 0,
                    "minH": 3,
                    "minW": 3
                  },
                  "sm": {
                    "h": 5,
                    "w": 4,
                    "x": 0,
                    "y": 5,
                    "minH": 3,
                    "minW": 3
                  },
                  "xl": {
                    "h": 6,
                    "w": 6,
                    "x": 8,
                    "y": 8,
                    "minH": 3,
                    "minW": 3
                  },
                  "xxl": {
                    "h": 6,
                    "w": 6,
                    "x": 8,
                    "y": 8,
                    "minH": 3,
                    "minW": 3
                  }
                },
                "quantity": "results.properties.optoelectronic.solar_cell.absorber_fabrication",
                "scale": "linear",
                "showinput": true
              },
              {
                "title": "Band gap",
                "type": "histogram",
                "layout": {
                  "lg": {
                    "h": 4,
                    "w": 12,
                    "x": 12,
                    "y": 16,
                    "minH": 3,
                    "minW": 8
                  },
                  "md": {
                    "h": 3,
                    "w": 8,
                    "x": 10,
                    "y": 14,
                    "minH": 3,
                    "minW": 8
                  },
                  "sm": {
                    "h": 3,
                    "w": 8,
                    "x": 4,
                    "y": 10,
                    "minH": 3,
                    "minW": 8
                  },
                  "xl": {
                    "h": 3,
                    "w": 8,
                    "x": 0,
                    "y": 8,
                    "minH": 3,
                    "minW": 8
                  },
                  "xxl": {
                    "h": 3,
                    "w": 8,
                    "x": 0,
                    "y": 11,
                    "minH": 3,
                    "minW": 8
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
                  "lg": {
                    "h": 6,
                    "w": 6,
                    "x": 18,
                    "y": 6,
                    "minH": 3,
                    "minW": 3
                  },
                  "md": {
                    "h": 6,
                    "w": 5,
                    "x": 0,
                    "y": 14,
                    "minH": 3,
                    "minW": 3
                  },
                  "sm": {
                    "h": 5,
                    "w": 4,
                    "x": 4,
                    "y": 5,
                    "minH": 3,
                    "minW": 3
                  },
                  "xl": {
                    "h": 6,
                    "w": 5,
                    "x": 25,
                    "y": 8,
                    "minH": 3,
                    "minW": 3
                  },
                  "xxl": {
                    "h": 6,
                    "w": 6,
                    "x": 20,
                    "y": 8,
                    "minH": 3,
                    "minW": 3
                  }
                },
                "quantity": "results.properties.optoelectronic.solar_cell.electron_transport_layer",
                "scale": "linear",
                "showinput": true
              },
              {
                "type": "terms",
                "layout": {
                  "lg": {
                    "h": 6,
                    "w": 6,
                    "x": 12,
                    "y": 6,
                    "minH": 3,
                    "minW": 3
                  },
                  "md": {
                    "h": 6,
                    "w": 5,
                    "x": 5,
                    "y": 14,
                    "minH": 3,
                    "minW": 3
                  },
                  "sm": {
                    "h": 5,
                    "w": 4,
                    "x": 8,
                    "y": 5,
                    "minH": 3,
                    "minW": 3
                  },
                  "xl": {
                    "h": 6,
                    "w": 5,
                    "x": 20,
                    "y": 8,
                    "minH": 3,
                    "minW": 3
                  },
                  "xxl": {
                    "h": 6,
                    "w": 6,
                    "x": 26,
                    "y": 8,
                    "minH": 3,
                    "minW": 3
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
          },
          "search_syntaxes": {
            "exclude": [
              "free_text"
            ]
          }
        },
        "mofs": {
          "label": "Metal-Organic Frameworks",
          "path": "mofs",
          "resource": "entries",
          "category": "Use Cases",
          "description": "Search metal-organic frameworks (MOFs)",
          "readme": "This page allows you to search **metal-organic framework (MOF) data** within NOMAD. The filter menu on the left and the shown default columns are specifically designed for MOF exploration. The dashboard directly shows useful interactive statistics about the data.",
          "pagination": {
            "order_by": "upload_create_time",
            "order": "desc",
            "page_size": 20
          },
          "columns": {
            "options": {
              "results.material.chemical_formula_iupac": {
                "label": "Formula",
                "align": "left"
              },
              "mainfile": {
                "align": "left"
              },
              "upload_create_time": {
                "label": "Upload time",
                "align": "left"
              },
              "authors": {
                "align": "left"
              },
              "comment": {
                "align": "left"
              },
              "datasets": {
                "align": "left"
              },
              "published": {
                "label": "Access",
                "align": "left"
              }
            },
            "selected": [
              "results.material.chemical_formula_iupac",
              "mainfile",
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
                "label": "Structure",
                "level": 1,
                "size": "s"
              },
              "electronic": {
                "label": "Electronic Properties",
                "level": 0,
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
                  "lg": {
                    "h": 9,
                    "w": 15,
                    "x": 0,
                    "y": 0,
                    "minH": 3,
                    "minW": 3
                  },
                  "md": {
                    "h": 8,
                    "w": 11,
                    "x": 0,
                    "y": 0,
                    "minH": 3,
                    "minW": 3
                  },
                  "sm": {
                    "h": 6,
                    "w": 9,
                    "x": 0,
                    "y": 0,
                    "minH": 3,
                    "minW": 3
                  },
                  "xl": {
                    "h": 9,
                    "w": 19,
                    "x": 0,
                    "y": 0,
                    "minH": 3,
                    "minW": 3
                  },
                  "xxl": {
                    "h": 10,
                    "w": 25,
                    "x": 0,
                    "y": 0,
                    "minH": 3,
                    "minW": 3
                  }
                },
                "quantity": "results.material.elements",
                "scale": "linear"
              },
              {
                "title": "SBU type",
                "type": "terms",
                "layout": {
                  "lg": {
                    "h": 9,
                    "w": 9,
                    "x": 15,
                    "y": 0,
                    "minH": 3,
                    "minW": 3
                  },
                  "md": {
                    "h": 8,
                    "w": 7,
                    "x": 11,
                    "y": 0,
                    "minH": 3,
                    "minW": 3
                  },
                  "sm": {
                    "h": 6,
                    "w": 3,
                    "x": 9,
                    "y": 0,
                    "minH": 3,
                    "minW": 3
                  },
                  "xl": {
                    "h": 9,
                    "w": 11,
                    "x": 19,
                    "y": 0,
                    "minH": 3,
                    "minW": 3
                  },
                  "xxl": {
                    "h": 10,
                    "w": 11,
                    "x": 25,
                    "y": 0,
                    "minH": 3,
                    "minW": 3
                  }
                },
                "quantity": "results.material.topology.sbu_type",
                "scale": "linear",
                "showinput": true
              },
              {
                "type": "histogram",
                "layout": {
                  "lg": {
                    "h": 5,
                    "w": 12,
                    "x": 0,
                    "y": 9,
                    "minH": 3,
                    "minW": 3
                  },
                  "md": {
                    "h": 4,
                    "w": 9,
                    "x": 0,
                    "y": 8,
                    "minH": 3,
                    "minW": 3
                  },
                  "sm": {
                    "h": 3,
                    "w": 6,
                    "x": 0,
                    "y": 6,
                    "minH": 3,
                    "minW": 3
                  },
                  "xl": {
                    "h": 5,
                    "w": 15,
                    "x": 0,
                    "y": 9,
                    "minH": 3,
                    "minW": 3
                  },
                  "xxl": {
                    "h": 6,
                    "w": 19,
                    "x": 0,
                    "y": 10,
                    "minH": 3,
                    "minW": 3
                  }
                },
                "quantity": "results.material.topology.pore_limiting_diameter",
                "scale": "linear",
                "autorange": false,
                "showinput": true,
                "nbins": 30
              },
              {
                "type": "histogram",
                "layout": {
                  "lg": {
                    "h": 5,
                    "w": 12,
                    "x": 0,
                    "y": 14,
                    "minH": 3,
                    "minW": 3
                  },
                  "md": {
                    "h": 4,
                    "w": 9,
                    "x": 9,
                    "y": 8,
                    "minH": 3,
                    "minW": 3
                  },
                  "sm": {
                    "h": 3,
                    "w": 6,
                    "x": 6,
                    "y": 6,
                    "minH": 3,
                    "minW": 3
                  },
                  "xl": {
                    "h": 5,
                    "w": 15,
                    "x": 0,
                    "y": 14,
                    "minH": 3,
                    "minW": 3
                  },
                  "xxl": {
                    "h": 6,
                    "w": 17,
                    "x": 19,
                    "y": 10,
                    "minH": 3,
                    "minW": 3
                  }
                },
                "quantity": "results.material.topology.largest_cavity_diameter",
                "scale": "linear",
                "autorange": false,
                "showinput": true,
                "nbins": 30
              },
              {
                "type": "histogram",
                "layout": {
                  "lg": {
                    "h": 5,
                    "w": 12,
                    "x": 11,
                    "y": 9,
                    "minH": 3,
                    "minW": 3
                  },
                  "md": {
                    "h": 4,
                    "w": 9,
                    "x": 0,
                    "y": 12,
                    "minH": 3,
                    "minW": 3
                  },
                  "sm": {
                    "h": 3,
                    "w": 6,
                    "x": 0,
                    "y": 9,
                    "minH": 3,
                    "minW": 3
                  },
                  "xl": {
                    "h": 5,
                    "w": 15,
                    "x": 15,
                    "y": 9,
                    "minH": 3,
                    "minW": 3
                  },
                  "xxl": {
                    "h": 6,
                    "w": 19,
                    "x": 0,
                    "y": 16,
                    "minH": 3,
                    "minW": 3
                  }
                },
                "quantity": "results.material.topology.accessible_surface_area",
                "scale": "linear",
                "autorange": false,
                "showinput": true,
                "nbins": 30
              },
              {
                "type": "histogram",
                "layout": {
                  "lg": {
                    "h": 5,
                    "w": 12,
                    "x": 11,
                    "y": 14,
                    "minH": 3,
                    "minW": 3
                  },
                  "md": {
                    "h": 4,
                    "w": 9,
                    "x": 9,
                    "y": 12,
                    "minH": 3,
                    "minW": 3
                  },
                  "sm": {
                    "h": 3,
                    "w": 6,
                    "x": 6,
                    "y": 9,
                    "minH": 3,
                    "minW": 3
                  },
                  "xl": {
                    "h": 5,
                    "w": 15,
                    "x": 15,
                    "y": 14,
                    "minH": 3,
                    "minW": 3
                  },
                  "xxl": {
                    "h": 6,
                    "w": 17,
                    "x": 19,
                    "y": 16,
                    "minH": 3,
                    "minW": 3
                  }
                },
                "quantity": "results.material.topology.void_fraction",
                "scale": "linear",
                "autorange": false,
                "showinput": true,
                "nbins": 30
              }
            ]
          },
          "filters_locked": {
            "results.material.topology.label": "MOF"
          },
          "search_syntaxes": {
            "exclude": [
              "free_text"
            ]
          }
        },
        "heterogeneouscatalyst": {
          "label": "Heterogeneous Catalysis",
          "path": "heterogeneouscatalyst",
          "resource": "entries",
          "category": "Use Cases",
          "description": "Search heterogeneous catalysts",
          "readme": "This page allows you to search **catalyst and catalysis data** within NOMAD. The filter menu on the left and the shown default columns are specifically designed for Heterogeneous Catalyst exploration. The dashboard directly shows useful interactive statistics about the data.",
          "pagination": {
            "order_by": "upload_create_time",
            "order": "desc",
            "page_size": 20
          },
          "columns": {
            "options": {
              "results.material.elements": {
                "align": "left"
              },
              "results.properties.catalytic.catalyst_synthesis.catalyst_type": {
                "align": "left"
              },
              "results.properties.catalytic.catalyst_synthesis.catalyst_name": {
                "align": "left"
              },
              "results.properties.catalytic.catalyst_synthesis.preparation_method": {
                "label": "Preparation",
                "align": "left"
              },
              "results.properties.catalytic.catalyst_characterization.surface_area": {
                "label": "Surface area (m^2/g)",
                "align": "left",
                "format": {
                  "decimals": 2,
                  "mode": "standard"
                }
              },
              "results.properties.catalytic.reaction.name": {
                "label": "Reaction name",
                "align": "left"
              },
              "results.properties.catalytic.reaction.type": {
                "label": "Reaction class",
                "align": "left"
              },
              "results.properties.catalytic.reaction.reactants.name": {
                "label": "Reactants",
                "align": "left"
              },
              "results.properties.catalytic.reaction.products.name": {
                "label": "Products",
                "align": "left"
              },
              "references": {
                "align": "left"
              },
              "results.material.chemical_formula_hill": {
                "label": "Formula",
                "align": "left"
              },
              "results.material.structural_type": {
                "align": "left"
              },
              "results.eln.lab_ids": {
                "align": "left"
              },
              "results.eln.sections": {
                "align": "left"
              },
              "results.eln.methods": {
                "align": "left"
              },
              "results.eln.tags": {
                "align": "left"
              },
              "results.eln.instruments": {
                "align": "left"
              },
              "entry_name": {
                "label": "Name",
                "align": "left"
              },
              "entry_type": {
                "align": "left"
              },
              "mainfile": {
                "align": "left"
              },
              "upload_create_time": {
                "label": "Upload time",
                "align": "left"
              },
              "authors": {
                "align": "left"
              },
              "comment": {
                "align": "left"
              },
              "datasets": {
                "align": "left"
              },
              "published": {
                "label": "Access",
                "align": "left"
              }
            },
            "selected": [
              "entry_name",
              "results.properties.catalytic.reaction.name",
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
                  "lg": {
                    "h": 10,
                    "w": 16,
                    "x": 0,
                    "y": 6,
                    "minH": 8,
                    "minW": 12
                  },
                  "md": {
                    "h": 8,
                    "w": 12,
                    "x": 0,
                    "y": 5,
                    "minH": 8,
                    "minW": 12
                  },
                  "sm": {
                    "h": 8,
                    "w": 12,
                    "x": 0,
                    "y": 4,
                    "minH": 8,
                    "minW": 12
                  },
                  "xl": {
                    "h": 10,
                    "w": 16,
                    "x": 0,
                    "y": 6,
                    "minH": 8,
                    "minW": 12
                  },
                  "xxl": {
                    "h": 10,
                    "w": 16,
                    "x": 0,
                    "y": 6,
                    "minH": 8,
                    "minW": 12
                  }
                },
                "quantity": "results.material.elements",
                "scale": "linear"
              },
              {
                "type": "terms",
                "layout": {
                  "lg": {
                    "h": 6,
                    "w": 8,
                    "x": 8,
                    "y": 0,
                    "minH": 3,
                    "minW": 3
                  },
                  "md": {
                    "h": 5,
                    "w": 6,
                    "x": 6,
                    "y": 0,
                    "minH": 3,
                    "minW": 3
                  },
                  "sm": {
                    "h": 4,
                    "w": 4,
                    "x": 4,
                    "y": 0,
                    "minH": 4,
                    "minW": 3
                  },
                  "xl": {
                    "h": 6,
                    "w": 8,
                    "x": 8,
                    "y": 0,
                    "minH": 3,
                    "minW": 3
                  },
                  "xxl": {
                    "h": 6,
                    "w": 8,
                    "x": 8,
                    "y": 0,
                    "minH": 3,
                    "minW": 3
                  }
                },
                "quantity": "results.properties.catalytic.reaction.reactants.name",
                "scale": "linear",
                "showinput": true
              },
              {
                "type": "terms",
                "layout": {
                  "lg": {
                    "h": 6,
                    "w": 8,
                    "x": 0,
                    "y": 0,
                    "minH": 3,
                    "minW": 3
                  },
                  "md": {
                    "h": 5,
                    "w": 6,
                    "x": 0,
                    "y": 0,
                    "minH": 3,
                    "minW": 3
                  },
                  "sm": {
                    "h": 4,
                    "w": 4,
                    "x": 0,
                    "y": 0,
                    "minH": 3,
                    "minW": 3
                  },
                  "xl": {
                    "h": 6,
                    "w": 8,
                    "x": 0,
                    "y": 0,
                    "minH": 3,
                    "minW": 3
                  },
                  "xxl": {
                    "h": 6,
                    "w": 8,
                    "x": 0,
                    "y": 0,
                    "minH": 3,
                    "minW": 3
                  }
                },
                "quantity": "results.properties.catalytic.reaction.name",
                "scale": "linear",
                "showinput": true
              },
              {
                "type": "terms",
                "layout": {
                  "lg": {
                    "h": 6,
                    "w": 8,
                    "x": 16,
                    "y": 0,
                    "minH": 3,
                    "minW": 3
                  },
                  "md": {
                    "h": 5,
                    "w": 6,
                    "x": 12,
                    "y": 0,
                    "minH": 3,
                    "minW": 3
                  },
                  "sm": {
                    "h": 4,
                    "w": 4,
                    "x": 8,
                    "y": 0,
                    "minH": 3,
                    "minW": 3
                  },
                  "xl": {
                    "h": 6,
                    "w": 8,
                    "x": 16,
                    "y": 0,
                    "minH": 3,
                    "minW": 3
                  },
                  "xxl": {
                    "h": 6,
                    "w": 8,
                    "x": 16,
                    "y": 0,
                    "minH": 3,
                    "minW": 3
                  }
                },
                "quantity": "results.properties.catalytic.reaction.products.name",
                "scale": "linear",
                "showinput": true
              },
              {
                "type": "terms",
                "layout": {
                  "lg": {
                    "h": 5,
                    "w": 8,
                    "x": 16,
                    "y": 6,
                    "minH": 3,
                    "minW": 3
                  },
                  "md": {
                    "h": 4,
                    "w": 6,
                    "x": 12,
                    "y": 5,
                    "minH": 3,
                    "minW": 3
                  },
                  "sm": {
                    "h": 3,
                    "w": 6,
                    "x": 0,
                    "y": 12,
                    "minH": 3,
                    "minW": 3
                  },
                  "xl": {
                    "h": 5,
                    "w": 8,
                    "x": 16,
                    "y": 6,
                    "minH": 3,
                    "minW": 3
                  },
                  "xxl": {
                    "h": 5,
                    "w": 8,
                    "x": 16,
                    "y": 6,
                    "minH": 3,
                    "minW": 3
                  }
                },
                "quantity": "results.properties.catalytic.catalyst_synthesis.preparation_method",
                "scale": "linear",
                "showinput": true
              },
              {
                "type": "terms",
                "layout": {
                  "lg": {
                    "h": 5,
                    "w": 8,
                    "x": 16,
                    "y": 11,
                    "minH": 3,
                    "minW": 3
                  },
                  "md": {
                    "h": 4,
                    "w": 6,
                    "x": 12,
                    "y": 9,
                    "minH": 3,
                    "minW": 3
                  },
                  "sm": {
                    "h": 3,
                    "w": 6,
                    "x": 6,
                    "y": 12,
                    "minH": 3,
                    "minW": 3
                  },
                  "xl": {
                    "h": 5,
                    "w": 8,
                    "x": 16,
                    "y": 11,
                    "minH": 3,
                    "minW": 3
                  },
                  "xxl": {
                    "h": 5,
                    "w": 8,
                    "x": 16,
                    "y": 11,
                    "minH": 3,
                    "minW": 3
                  }
                },
                "quantity": "results.properties.catalytic.catalyst_synthesis.catalyst_type",
                "scale": "linear",
                "showinput": true
              },
              {
                "type": "histogram",
                "layout": {
                  "lg": {
                    "h": 5,
                    "w": 12,
                    "x": 12,
                    "y": 21,
                    "minH": 3,
                    "minW": 8
                  },
                  "md": {
                    "h": 3,
                    "w": 9,
                    "x": 9,
                    "y": 16,
                    "minH": 3,
                    "minW": 8
                  },
                  "sm": {
                    "h": 3,
                    "w": 6,
                    "x": 6,
                    "y": 18,
                    "minH": 3,
                    "minW": 6
                  },
                  "xl": {
                    "h": 4,
                    "w": 12,
                    "x": 12,
                    "y": 20,
                    "minH": 3,
                    "minW": 8
                  },
                  "xxl": {
                    "h": 4,
                    "w": 12,
                    "x": 12,
                    "y": 20,
                    "minH": 3,
                    "minW": 8
                  }
                },
                "quantity": "results.properties.catalytic.reaction.weight_hourly_space_velocity",
                "scale": "linear",
                "autorange": false,
                "showinput": false,
                "nbins": 30
              },
              {
                "type": "scatterplot",
                "layout": {
                  "lg": {
                    "h": 10,
                    "w": 12,
                    "x": 0,
                    "y": 16,
                    "minH": 3,
                    "minW": 8
                  },
                  "md": {
                    "h": 6,
                    "w": 9,
                    "x": 0,
                    "y": 13,
                    "minH": 3,
                    "minW": 8
                  },
                  "sm": {
                    "h": 6,
                    "w": 6,
                    "x": 0,
                    "y": 15,
                    "minH": 3,
                    "minW": 6
                  },
                  "xl": {
                    "h": 8,
                    "w": 12,
                    "x": 0,
                    "y": 16,
                    "minH": 3,
                    "minW": 8
                  },
                  "xxl": {
                    "h": 8,
                    "w": 12,
                    "x": 0,
                    "y": 16,
                    "minH": 6,
                    "minW": 8
                  }
                },
                "x": {
                  "title": "gas concentration (%)",
                  "quantity": "results.properties.catalytic.reaction.reactants[*].gas_concentration_in"
                },
                "y": {
                  "quantity": "results.properties.catalytic.reaction.temperature"
                },
                "markers": {
                  "color": {
                    "quantity": "results.properties.catalytic.reaction.reactants[*].name"
                  }
                },
                "size": 1000,
                "autorange": true
              },
              {
                "type": "histogram",
                "layout": {
                  "lg": {
                    "h": 5,
                    "w": 12,
                    "x": 12,
                    "y": 16,
                    "minH": 3,
                    "minW": 8
                  },
                  "md": {
                    "h": 3,
                    "w": 9,
                    "x": 9,
                    "y": 13,
                    "minH": 3,
                    "minW": 8
                  },
                  "sm": {
                    "h": 3,
                    "w": 6,
                    "x": 6,
                    "y": 15,
                    "minH": 3,
                    "minW": 6
                  },
                  "xl": {
                    "h": 4,
                    "w": 12,
                    "x": 12,
                    "y": 16,
                    "minH": 3,
                    "minW": 8
                  },
                  "xxl": {
                    "h": 4,
                    "w": 12,
                    "x": 12,
                    "y": 16,
                    "minH": 3,
                    "minW": 8
                  }
                },
                "quantity": "results.properties.catalytic.reaction.pressure",
                "scale": "linear",
                "autorange": false,
                "showinput": false,
                "nbins": 30
              },
              {
                "type": "scatterplot",
                "layout": {
                  "lg": {
                    "h": 8,
                    "w": 12,
                    "x": 0,
                    "y": 26,
                    "minH": 3,
                    "minW": 3
                  },
                  "md": {
                    "h": 6,
                    "w": 9,
                    "x": 0,
                    "y": 19,
                    "minH": 3,
                    "minW": 3
                  },
                  "sm": {
                    "h": 5,
                    "w": 6,
                    "x": 0,
                    "y": 21,
                    "minH": 3,
                    "minW": 3
                  },
                  "xl": {
                    "h": 9,
                    "w": 12,
                    "x": 0,
                    "y": 24,
                    "minH": 3,
                    "minW": 3
                  },
                  "xxl": {
                    "h": 9,
                    "w": 12,
                    "x": 0,
                    "y": 24,
                    "minH": 3,
                    "minW": 3
                  }
                },
                "x": {
                  "quantity": "results.properties.catalytic.reaction.temperature"
                },
                "y": {
                  "title": "Conversion (%)",
                  "quantity": "results.properties.catalytic.reaction.reactants[*].conversion"
                },
                "markers": {
                  "color": {
                    "quantity": "results.properties.catalytic.reaction.reactants[*].name"
                  }
                },
                "size": 1000,
                "autorange": true
              },
              {
                "type": "scatterplot",
                "layout": {
                  "lg": {
                    "h": 8,
                    "w": 12,
                    "x": 12,
                    "y": 26,
                    "minH": 3,
                    "minW": 3
                  },
                  "md": {
                    "h": 6,
                    "w": 9,
                    "x": 9,
                    "y": 19,
                    "minH": 3,
                    "minW": 3
                  },
                  "sm": {
                    "h": 5,
                    "w": 6,
                    "x": 6,
                    "y": 21,
                    "minH": 3,
                    "minW": 3
                  },
                  "xl": {
                    "h": 9,
                    "w": 12,
                    "x": 12,
                    "y": 24,
                    "minH": 3,
                    "minW": 3
                  },
                  "xxl": {
                    "h": 9,
                    "w": 12,
                    "x": 12,
                    "y": 24,
                    "minH": 3,
                    "minW": 3
                  }
                },
                "x": {
                  "quantity": "results.properties.catalytic.reaction.temperature"
                },
                "y": {
                  "title": "Selectivity (%)",
                  "quantity": "results.properties.catalytic.reaction.products[*].selectivity"
                },
                "markers": {
                  "color": {
                    "quantity": "results.properties.catalytic.reaction.products[*].name"
                  }
                },
                "size": 1000,
                "autorange": true
              },
              {
                "type": "scatterplot",
                "layout": {
                  "lg": {
                    "h": 8,
                    "w": 12,
                    "x": 0,
                    "y": 34,
                    "minH": 3,
                    "minW": 3
                  },
                  "md": {
                    "h": 6,
                    "w": 9,
                    "x": 0,
                    "y": 25,
                    "minH": 3,
                    "minW": 3
                  },
                  "sm": {
                    "h": 5,
                    "w": 6,
                    "x": 0,
                    "y": 26,
                    "minH": 3,
                    "minW": 3
                  },
                  "xl": {
                    "h": 9,
                    "w": 12,
                    "x": 0,
                    "y": 33,
                    "minH": 3,
                    "minW": 3
                  },
                  "xxl": {
                    "h": 9,
                    "w": 12,
                    "x": 0,
                    "y": 33,
                    "minH": 3,
                    "minW": 3
                  }
                },
                "x": {
                  "title": "Oxygen Conversion (%)",
                  "quantity": "results.properties.catalytic.reaction.reactants[? name=='molecular oxygen'].conversion"
                },
                "y": {
                  "title": "Acetic Acid Selectivity (%)",
                  "quantity": "results.properties.catalytic.reaction.products[? name=='acetic acid'].selectivity"
                },
                "markers": {
                  "color": {
                    "quantity": "results.properties.catalytic.reaction.name"
                  }
                },
                "size": 1000,
                "autorange": true
              },
              {
                "type": "scatterplot",
                "layout": {
                  "lg": {
                    "h": 8,
                    "w": 12,
                    "x": 12,
                    "y": 34,
                    "minH": 3,
                    "minW": 3
                  },
                  "md": {
                    "h": 6,
                    "w": 9,
                    "x": 9,
                    "y": 25,
                    "minH": 3,
                    "minW": 3
                  },
                  "sm": {
                    "h": 5,
                    "w": 6,
                    "x": 6,
                    "y": 26,
                    "minH": 3,
                    "minW": 3
                  },
                  "xl": {
                    "h": 9,
                    "w": 12,
                    "x": 12,
                    "y": 33,
                    "minH": 3,
                    "minW": 3
                  },
                  "xxl": {
                    "h": 9,
                    "w": 12,
                    "x": 12,
                    "y": 33,
                    "minH": 3,
                    "minW": 3
                  }
                },
                "x": {
                  "title": "Carbon Monoxide Conversion (%)",
                  "quantity": "results.properties.catalytic.reaction.reactants[? name=='carbon monoxide'].conversion"
                },
                "y": {
                  "title": "Ethanol Selectivity (%)",
                  "quantity": "results.properties.catalytic.reaction.products[? name=='ethanol'].selectivity"
                },
                "markers": {
                  "color": {
                    "quantity": "results.properties.catalytic.catalyst_synthesis.preparation_method"
                  }
                },
                "size": 1000,
                "autorange": true
              },
              {
                "type": "histogram",
                "layout": {
                  "lg": {
                    "h": 4,
                    "w": 12,
                    "x": 0,
                    "y": 42,
                    "minH": 3,
                    "minW": 8
                  },
                  "md": {
                    "h": 4,
                    "w": 9,
                    "x": 0,
                    "y": 31,
                    "minH": 3,
                    "minW": 8
                  },
                  "sm": {
                    "h": 4,
                    "w": 8,
                    "x": 0,
                    "y": 31,
                    "minH": 3,
                    "minW": 8
                  },
                  "xl": {
                    "h": 4,
                    "w": 12,
                    "x": 0,
                    "y": 42,
                    "minH": 3,
                    "minW": 8
                  },
                  "xxl": {
                    "h": 4,
                    "w": 12,
                    "x": 0,
                    "y": 42,
                    "minH": 3,
                    "minW": 8
                  }
                },
                "quantity": "results.properties.catalytic.catalyst_characterization.surface_area",
                "scale": "1/4",
                "autorange": false,
                "showinput": false,
                "nbins": 30
              },
              {
                "type": "terms",
                "layout": {
                  "lg": {
                    "h": 4,
                    "w": 12,
                    "x": 12,
                    "y": 42,
                    "minH": 3,
                    "minW": 3
                  },
                  "md": {
                    "h": 4,
                    "w": 9,
                    "x": 9,
                    "y": 31,
                    "minH": 3,
                    "minW": 3
                  },
                  "sm": {
                    "h": 4,
                    "w": 8,
                    "x": 8,
                    "y": 31,
                    "minH": 4,
                    "minW": 3
                  },
                  "xl": {
                    "h": 4,
                    "w": 12,
                    "x": 12,
                    "y": 42,
                    "minH": 3,
                    "minW": 3
                  },
                  "xxl": {
                    "h": 4,
                    "w": 12,
                    "x": 12,
                    "y": 42,
                    "minH": 3,
                    "minW": 3
                  }
                },
                "quantity": "results.properties.catalytic.catalyst_characterization.method",
                "scale": "linear",
                "showinput": true
              }
            ]
          },
          "filters_locked": {
            "quantities": "results.properties.catalytic"
          },
          "search_syntaxes": {
            "exclude": [
              "free_text"
            ]
          }
        }
      }
    },
    "north": {
      "enabled": true
    },
    "example_uploads": {}
  },
  "plugins": {
    "entry_points": {
      "include": [
        "schema/simulation/run",
        "schema/simulation/workflow",
        "parsers/vasp"
      ],
      "exclude": null,
      "options": {
        "normalizers/simulation/band_structure": {
          "plugin_type": "normalizer",
          "id": "normalizers/simulation/band_structure",
          "name": "bandstructurenormalizer",
          "description": "This is the normalizer for band structure in NOMAD.\n"
        },
        "normalizers/simulation/dos": {
          "plugin_type": "normalizer",
          "id": "normalizers/simulation/dos",
          "name": "dosnormalizer",
          "description": "This is the normalizer for DOS in NOMAD.\n"
        },
        "normalizers/simulation/soap": {
          "plugin_type": "normalizer",
          "id": "normalizers/simulation/soap",
          "name": "soapnormalizer",
          "description": "This is the normalizer for SOAP in NOMAD.\n"
        },
        "normalizers/simulation/spectra": {
          "plugin_type": "normalizer",
          "id": "normalizers/simulation/spectra",
          "name": "spectranormalizer",
          "description": "This is the normalizer for spectra in NOMAD.\n"
        },
        "normalizers/simulation/system": {
          "plugin_type": "normalizer",
          "id": "normalizers/simulation/system",
          "name": "systemnormalizer",
          "description": "This is the normalizer for system in NOMAD.\n"
        },
        "normalizers/simulation/workflow": {
          "plugin_type": "normalizer",
          "id": "normalizers/simulation/workflow",
          "name": "simulationworkflownormalizer",
          "description": "This is the normalizer for simulation workflows in NOMAD.\n"
        },
        "parsers/abacus": {
          "plugin_type": "parser",
          "id": "parsers/abacus",
          "name": "parsers/abacus",
          "plugin_source_code_url": "https://github.com/nomad-coe/electronic-parsers/tree/develop/electronicparsers/abacus"
        },
        "parsers/abinit": {
          "plugin_type": "parser",
          "id": "parsers/abinit",
          "name": "parsers/abinit",
          "plugin_source_code_url": "https://github.com/nomad-coe/electronic-parsers/tree/develop/electronicparsers/abinit"
        },
        "parsers/aflow": {
          "plugin_type": "parser",
          "id": "parsers/aflow",
          "name": "parsers/aflow",
          "plugin_source_code_url": "https://github.com/nomad-coe/workflow-parsers/tree/master/workflowparsers/aflow"
        },
        "parsers/amber": {
          "plugin_type": "parser",
          "id": "parsers/amber",
          "name": "parsers/amber",
          "plugin_source_code_url": "https://github.com/nomad-coe/atomistic-parsers/tree/develop/atomisticparsers/amber"
        },
        "parsers/ams": {
          "plugin_type": "parser",
          "id": "parsers/ams",
          "name": "parsers/ams",
          "plugin_source_code_url": "https://github.com/nomad-coe/electronic-parsers/tree/develop/electronicparsers/ams"
        },
        "parsers/asap": {
          "plugin_type": "parser",
          "id": "parsers/asap",
          "name": "parsers/asap",
          "plugin_source_code_url": "https://github.com/nomad-coe/atomistic-parsers/tree/develop/atomisticparsers/asap"
        },
        "parsers/asr": {
          "plugin_type": "parser",
          "id": "parsers/asr",
          "name": "parsers/asr",
          "plugin_source_code_url": "https://github.com/nomad-coe/workflow-parsers/tree/master/workflowparsers/asr"
        },
        "parsers/atk": {
          "plugin_type": "parser",
          "id": "parsers/atk",
          "name": "parsers/atk",
          "plugin_source_code_url": "https://github.com/nomad-coe/electronic-parsers/tree/develop/electronicparsers/atk"
        },
        "parsers/atomate": {
          "plugin_type": "parser",
          "id": "parsers/atomate",
          "name": "parsers/atomate",
          "plugin_source_code_url": "https://github.com/nomad-coe/workflow-parsers/tree/master/workflowparsers/automate"
        },
        "parsers/bigdft": {
          "plugin_type": "parser",
          "id": "parsers/bigdft",
          "name": "parsers/bigdft",
          "plugin_source_code_url": "https://github.com/nomad-coe/electronic-parsers/tree/develop/electronicparsers/bigdft"
        },
        "parsers/bopfox": {
          "plugin_type": "parser",
          "id": "parsers/bopfox",
          "name": "parsers/bopfox",
          "plugin_source_code_url": "https://github.com/nomad-coe/atomistic-parsers/tree/develop/atomisticparsers/bobfox"
        },
        "parsers/castep": {
          "plugin_type": "parser",
          "id": "parsers/castep",
          "name": "parsers/castep",
          "plugin_source_code_url": "https://github.com/nomad-coe/electronic-parsers/tree/develop/electronicparsers/castep"
        },
        "parsers/charmm": {
          "plugin_type": "parser",
          "id": "parsers/charmm",
          "name": "parsers/charmm",
          "plugin_source_code_url": "https://github.com/nomad-coe/electronic-parsers/tree/develop/electronicparsers/charmm"
        },
        "parsers/chemotion/chemotion": {
          "plugin_type": "parser",
          "id": "parsers/chemotion/chemotion",
          "name": "parsers/chemotion"
        },
        "parsers/cp2k": {
          "plugin_type": "parser",
          "id": "parsers/cp2k",
          "name": "parsers/cp2k",
          "plugin_source_code_url": "https://github.com/nomad-coe/electronic-parsers/tree/develop/electronicparsers/cp2k"
        },
        "parsers/cpmd": {
          "plugin_type": "parser",
          "id": "parsers/cpmd",
          "name": "parsers/cpmd",
          "plugin_source_code_url": "https://github.com/nomad-coe/electronic-parsers/tree/develop/electronicparsers/cpmd"
        },
        "parsers/crystal": {
          "plugin_type": "parser",
          "id": "parsers/crystal",
          "name": "parsers/crystal",
          "plugin_source_code_url": "https://github.com/nomad-coe/electronic-parsers/tree/develop/electronicparsers/crystal"
        },
        "parsers/dftbplus": {
          "plugin_type": "parser",
          "id": "parsers/dftbplus",
          "name": "parsers/dftbplus",
          "plugin_source_code_url": "https://github.com/nomad-coe/atomistic-parsers/tree/develop/atomisticparsers/dftplus"
        },
        "parsers/dlpoly": {
          "plugin_type": "parser",
          "id": "parsers/dlpoly",
          "name": "parsers/dl-poly",
          "plugin_source_code_url": "https://github.com/nomad-coe/atomistic-parsers/tree/develop/atomisticparsers/dlpoly"
        },
        "parsers/dmol3": {
          "plugin_type": "parser",
          "id": "parsers/dmol3",
          "name": "parsers/dmol",
          "plugin_source_code_url": "https://github.com/nomad-coe/electronic-parsers/tree/develop/electronicparsers/dmol3"
        },
        "parsers/edmft": {
          "plugin_type": "parser",
          "id": "parsers/edmft",
          "name": "parsers/edmft",
          "plugin_source_code_url": "https://github.com/nomad-coe/electronic-parsers/tree/develop/electronicparsers/edmft"
        },
        "parsers/eelsdbparser": {
          "plugin_type": "parser",
          "id": "parsers/eelsdbparser",
          "name": "parsers/eels",
          "plugin_source_code_url": "https://github.com/nomad-coe/nomad-parser-eelsdb"
        },
        "parsers/elabftw/elabftw": {
          "plugin_type": "parser",
          "id": "parsers/elabftw/elabftw",
          "name": "parsers/elabftw",
          "plugin_source_code_url": "https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-FAIR/-/tree/develop/nomad/parsing/elabftw"
        },
        "parsers/elastic": {
          "plugin_type": "parser",
          "id": "parsers/elastic",
          "name": "parsers/elastic",
          "plugin_source_code_url": "https://github.com/nomad-coe/workflow-parsers/tree/master/workflowparsers/elastic"
        },
        "parsers/elk": {
          "plugin_type": "parser",
          "id": "parsers/elk",
          "name": "parsers/elk",
          "plugin_source_code_url": "https://github.com/nomad-coe/electronic-parsers/tree/develop/electronicparsers/elk"
        },
        "parsers/exciting": {
          "plugin_type": "parser",
          "id": "parsers/exciting",
          "name": "parsers/exciting",
          "plugin_source_code_url": "https://github.com/nomad-coe/electronic-parsers/tree/develop/electronicparsers/exciting"
        },
        "parsers/fhi-aims": {
          "plugin_type": "parser",
          "id": "parsers/fhi-aims",
          "name": "parsers/fhi-aims",
          "plugin_source_code_url": "https://github.com/nomad-coe/electronic-parsers/tree/develop/electronicparsers/fhiaims"
        },
        "parsers/fhivibes": {
          "plugin_type": "parser",
          "id": "parsers/fhivibes",
          "name": "parsers/fhi-vibes",
          "plugin_source_code_url": "https://github.com/nomad-coe/workflow-parsers/tree/master/workflowparsers/fhivibes"
        },
        "parsers/fleur": {
          "plugin_type": "parser",
          "id": "parsers/fleur",
          "name": "parsers/fleur",
          "plugin_source_code_url": "https://github.com/nomad-coe/electronic-parsers/tree/develop/electronicparsers/fleur"
        },
        "parsers/fplo": {
          "plugin_type": "parser",
          "id": "parsers/fplo",
          "name": "parsers/fplo",
          "plugin_source_code_url": "https://github.com/nomad-coe/electronic-parsers/tree/develop/electronicparsers/fplo"
        },
        "parsers/gamess": {
          "plugin_type": "parser",
          "id": "parsers/gamess",
          "name": "parsers/gamess",
          "plugin_source_code_url": "https://github.com/nomad-coe/electronic-parsers/tree/develop/electronicparsers/gamess"
        },
        "parsers/gaussian": {
          "plugin_type": "parser",
          "id": "parsers/gaussian",
          "name": "parsers/gaussian",
          "plugin_source_code_url": "https://github.com/nomad-coe/electronic-parsers/tree/develop/electronicparsers/gaussian"
        },
        "parsers/gpaw": {
          "plugin_type": "parser",
          "id": "parsers/gpaw",
          "name": "parsers/gpaw",
          "plugin_source_code_url": "https://github.com/nomad-coe/electronic-parsers/tree/develop/electronicparsers/gpaw"
        },
        "parsers/gromacs": {
          "plugin_type": "parser",
          "id": "parsers/gromacs",
          "name": "parsers/gromacs",
          "plugin_source_code_url": "https://github.com/nomad-coe/atomistic-parsers/tree/develop/atomisticparsers/gromacs"
        },
        "parsers/gromos": {
          "plugin_type": "parser",
          "id": "parsers/gromos",
          "name": "parsers/gromos",
          "plugin_source_code_url": "https://github.com/nomad-coe/atomistic-parsers/tree/develop/atomisticparsers/gromos"
        },
        "parsers/gulp": {
          "plugin_type": "parser",
          "id": "parsers/gulp",
          "name": "parsers/gulp",
          "plugin_source_code_url": "https://github.com/nomad-coe/atomistic-parsers/tree/develop/atomisticparsers/gulp"
        },
        "parsers/h5md": {
          "plugin_type": "parser",
          "id": "parsers/h5md",
          "name": "parsers/h5md"
        },
        "parsers/lammps": {
          "plugin_type": "parser",
          "id": "parsers/lammps",
          "name": "parsers/lammps",
          "plugin_source_code_url": "https://github.com/nomad-coe/atomistic-parsers/tree/develop/atomisticparsers/lammps"
        },
        "parsers/libatoms": {
          "plugin_type": "parser",
          "id": "parsers/libatoms",
          "name": "parsers/lib-atoms",
          "plugin_source_code_url": "https://github.com/nomad-coe/atomistic-parsers/tree/develop/atomisticparsers/libatoms"
        },
        "parsers/lobster": {
          "plugin_type": "parser",
          "id": "parsers/lobster",
          "name": "parsers/lobster",
          "plugin_source_code_url": "https://github.com/nomad-coe/workflow-parsers/tree/master/workflowparsers/lobster"
        },
        "parsers/magres": {
          "plugin_type": "parser",
          "id": "parsers/magres",
          "name": "parsers/magres",
          "plugin_source_code_url": "https://github.com/nomad-coe/electronic-parsers/tree/develop/electronicparsers/magres"
        },
        "parsers/molcas": {
          "plugin_type": "parser",
          "id": "parsers/molcas",
          "name": "parsers/molcas",
          "plugin_source_code_url": "https://github.com/nomad-coe/electronic-parsers/tree/develop/electronicparsers/molcas"
        },
        "parsers/mopac": {
          "plugin_type": "parser",
          "id": "parsers/mopac",
          "name": "parsers/mopac",
          "plugin_source_code_url": "https://github.com/nomad-coe/electronic-parsers/tree/develop/electronicparsers/mopac"
        },
        "parsers/namd": {
          "plugin_type": "parser",
          "id": "parsers/namd",
          "name": "parsers/namd",
          "plugin_source_code_url": "https://github.com/nomad-coe/atomistic-parsers/tree/develop/atomisticparsers/namd"
        },
        "parsers/nexus": {
          "plugin_type": "parser",
          "id": "parsers/nexus",
          "name": "parsers/nexus",
          "plugin_source_code_url": "https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-FAIR/-/tree/develop/nomad/parsing/nexus"
        },
        "parsers/nwchem": {
          "plugin_type": "parser",
          "id": "parsers/nwchem",
          "name": "parsers/nwchem",
          "plugin_source_code_url": "https://github.com/nomad-coe/electronic-parsers/tree/develop/electronicparsers/nwchem"
        },
        "parsers/ocean": {
          "plugin_type": "parser",
          "id": "parsers/ocean",
          "name": "parsers/ocean",
          "plugin_source_code_url": "https://github.com/nomad-coe/electronic-parsers/tree/develop/electronicparsers/ocean"
        },
        "parsers/octopus": {
          "plugin_type": "parser",
          "id": "parsers/octopus",
          "name": "parsers/octopus",
          "plugin_source_code_url": "https://github.com/nomad-coe/electronic-parsers/tree/develop/electronicparsers/octopus"
        },
        "parsers/onetep": {
          "plugin_type": "parser",
          "id": "parsers/onetep",
          "name": "parsers/onetep",
          "plugin_source_code_url": "https://github.com/nomad-coe/electronic-parsers/tree/develop/electronicparsers/onetep"
        },
        "parsers/openkim": {
          "plugin_type": "parser",
          "id": "parsers/openkim",
          "name": "parsers/openkim",
          "plugin_source_code_url": "https://github.com/nomad-coe/electronic-parsers/tree/develop/databaseparsers/openkim"
        },
        "parsers/openmx": {
          "plugin_type": "parser",
          "id": "parsers/openmx",
          "name": "parsers/openmx",
          "plugin_source_code_url": "https://github.com/nomad-coe/electronic-parsers/tree/develop/electronicparsers/openmx"
        },
        "parsers/orca": {
          "plugin_type": "parser",
          "id": "parsers/orca",
          "name": "parsers/orca",
          "plugin_source_code_url": "https://github.com/nomad-coe/electronic-parsers/tree/develop/electronicparsers/orca"
        },
        "parsers/phonopy": {
          "plugin_type": "parser",
          "id": "parsers/phonopy",
          "name": "parsers/phonopy",
          "plugin_source_code_url": "https://github.com/nomad-coe/workflow-parsers/tree/master/workflowparsers/phonopy"
        },
        "parsers/psi4": {
          "plugin_type": "parser",
          "id": "parsers/psi4",
          "name": "parsers/psi4",
          "plugin_source_code_url": "https://github.com/nomad-coe/electronic-parsers/tree/develop/electronicparsers/psi4"
        },
        "parsers/qball": {
          "plugin_type": "parser",
          "id": "parsers/qball",
          "name": "parsers/qball",
          "plugin_source_code_url": "https://github.com/nomad-coe/electronic-parsers/tree/develop/electronicparsers/qball"
        },
        "parsers/qbox": {
          "plugin_type": "parser",
          "id": "parsers/qbox",
          "name": "parsers/qbox",
          "plugin_source_code_url": "https://github.com/nomad-coe/electronic-parsers/tree/develop/electronicparsers/qbox"
        },
        "parsers/quantum_espresso_epw": {
          "plugin_type": "parser",
          "id": "parsers/quantum_espresso_epw",
          "name": "parsers/quantumespressoepw",
          "plugin_source_code_url": "https://github.com/nomad-coe/workflow-parsers/tree/master/workflowparsers/quantum_espresso_epw"
        },
        "parsers/quantum_espresso_phonon": {
          "plugin_type": "parser",
          "id": "parsers/quantum_espresso_phonon",
          "name": "parsers/quantumespressophonon",
          "plugin_source_code_url": "https://github.com/nomad-coe/workflow-parsers/tree/master/workflowparsers/quantum_espresso_phonon"
        },
        "parsers/quantum_espresso_xspectra": {
          "plugin_type": "parser",
          "id": "parsers/quantum_espresso_xspectra",
          "name": "parsers/quantumespressoxspectra",
          "plugin_source_code_url": "https://github.com/nomad-coe/workflow-parsers/tree/master/workflowparsers/quantum_espresso_xpectra"
        },
        "parsers/quantumespresso": {
          "plugin_type": "parser",
          "id": "parsers/quantumespresso",
          "name": "parsers/quantumespresso",
          "plugin_source_code_url": "https://github.com/nomad-coe/electronic-parsers/tree/develop/electronicparsers/quantumespresso"
        },
        "parsers/siesta": {
          "plugin_type": "parser",
          "id": "parsers/siesta",
          "name": "parsers/siesta",
          "plugin_source_code_url": "https://github.com/nomad-coe/electronic-parsers/tree/develop/electronicparsers/siesta"
        },
        "parsers/soliddmft": {
          "plugin_type": "parser",
          "id": "parsers/soliddmft",
          "name": "parsers/soliddmft",
          "plugin_source_code_url": "https://github.com/nomad-coe/electronic-parsers/tree/develop/electronicparsers/soliddmft"
        },
        "parsers/tbstudio": {
          "plugin_type": "parser",
          "id": "parsers/tbstudio",
          "name": "parsers/tbstudio"
        },
        "parsers/tinker": {
          "plugin_type": "parser",
          "id": "parsers/tinker",
          "name": "parsers/tinker",
          "plugin_source_code_url": "https://github.com/nomad-coe/atomistic-parsers/tree/develop/atomisticparsers/tinker"
        },
        "parsers/turbomole": {
          "plugin_type": "parser",
          "id": "parsers/turbomole",
          "name": "parsers/turbomole",
          "plugin_source_code_url": "https://github.com/nomad-coe/electronic-parsers/tree/develop/electronicparsers/turbomole"
        },
        "parsers/vasp": {
          "plugin_type": "parser",
          "id": "parsers/vasp",
          "name": "parsers/vasp",
          "plugin_source_code_url": "https://github.com/nomad-coe/electronic-parsers/tree/develop/electronicparsers/vasp"
        },
        "parsers/w2dynamics": {
          "plugin_type": "parser",
          "id": "parsers/w2dynamics",
          "name": "parsers/w2dynamics",
          "plugin_source_code_url": "https://github.com/nomad-coe/electronic-parsers/tree/develop/electronicparsers/w2dynamics"
        },
        "parsers/wannier90": {
          "plugin_type": "parser",
          "id": "parsers/wannier90",
          "name": "parsers/wannier90",
          "plugin_source_code_url": "https://github.com/nomad-coe/electronic-parsers/tree/develop/electronicparsers/wannier90"
        },
        "parsers/wien2k": {
          "plugin_type": "parser",
          "id": "parsers/wien2k",
          "name": "parsers/wien2k",
          "plugin_source_code_url": "https://github.com/nomad-coe/electronic-parsers/tree/develop/electronicparsers/wien2k"
        },
        "parsers/xtb": {
          "plugin_type": "parser",
          "id": "parsers/xtb",
          "name": "parsers/xtb",
          "plugin_source_code_url": "https://github.com/nomad-coe/atomistic-parsers/tree/develop/atomisticparsers/xtb"
        },
        "parsers/yambo": {
          "plugin_type": "parser",
          "id": "parsers/yambo",
          "name": "parsers/yambo",
          "plugin_source_code_url": "https://github.com/nomad-coe/electronic-parsers/tree/develop/electronicparsers/yambo"
        },
        "schema/nomad-perovskite-solar-cells-database/perovskite_solar_cell_database": {
          "plugin_type": "schema",
          "id": "schema/nomad-perovskite-solar-cells-database/perovskite_solar_cell_database",
          "name": "perovskite_solar_cell_database",
          "description": "A NOMAD plugin containing the schema for the Perovskite Solar Cell Database."
        },
        "schema/simulation/run": {
          "plugin_type": "schema",
          "id": "schema/simulation/run",
          "name": "runschema",
          "description": "Run schema plugin for NOMAD.\n"
        },
        "schema/simulation/workflow": {
          "plugin_type": "schema",
          "id": "schema/simulation/workflow",
          "name": "simulationworkflowschema",
          "description": "This is a collection of schemas for various types of simulation workflows.\n"
        }
      }
    },
    "plugin_packages": {}
  }
}
