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
                "label": "Entry type",
                "align": "left"
              },
              "entry_create_time": {
                "label": "Entry creation time",
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
                "label": "Min energy",
                "align": "left"
              },
              "results.properties.spectroscopic.spectra.provenance.eels.max_energy": {
                "label": "Max energy",
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
                  "lg": {
                    "h": 8,
                    "w": 12,
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
                    "y": 5,
                    "minH": 8,
                    "minW": 12
                  },
                  "xl": {
                    "h": 8,
                    "w": 12,
                    "x": 0,
                    "y": 5,
                    "minH": 8,
                    "minW": 12
                  },
                  "xxl": {
                    "h": 8,
                    "w": 12,
                    "x": 0,
                    "y": 5,
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
                    "w": 6,
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
                    "h": 5,
                    "w": 4,
                    "x": 0,
                    "y": 0,
                    "minH": 3,
                    "minW": 3
                  },
                  "xl": {
                    "h": 5,
                    "w": 6,
                    "x": 0,
                    "y": 0,
                    "minH": 3,
                    "minW": 3
                  },
                  "xxl": {
                    "h": 5,
                    "w": 6,
                    "x": 6,
                    "y": 0,
                    "minH": 3,
                    "minW": 3
                  }
                },
                "quantity": "results.properties.catalytic.reactivity.reactants.name",
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
                    "h": 5,
                    "w": 4,
                    "x": 8,
                    "y": 0,
                    "minH": 3,
                    "minW": 3
                  },
                  "xl": {
                    "h": 5,
                    "w": 6,
                    "x": 12,
                    "y": 0,
                    "minH": 3,
                    "minW": 3
                  },
                  "xxl": {
                    "h": 5,
                    "w": 6,
                    "x": 0,
                    "y": 0,
                    "minH": 3,
                    "minW": 3
                  }
                },
                "quantity": "results.properties.catalytic.reactivity.reaction_name",
                "scale": "linear",
                "showinput": true
              },
              {
                "type": "terms",
                "layout": {
                  "lg": {
                    "h": 6,
                    "w": 6,
                    "x": 6,
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
                    "h": 5,
                    "w": 4,
                    "x": 4,
                    "y": 0,
                    "minH": 3,
                    "minW": 3
                  },
                  "xl": {
                    "h": 5,
                    "w": 6,
                    "x": 6,
                    "y": 0,
                    "minH": 3,
                    "minW": 3
                  },
                  "xxl": {
                    "h": 5,
                    "w": 6,
                    "x": 12,
                    "y": 0,
                    "minH": 3,
                    "minW": 3
                  }
                },
                "quantity": "results.properties.catalytic.reactivity.products.name",
                "scale": "linear",
                "showinput": true
              },
              {
                "type": "terms",
                "layout": {
                  "lg": {
                    "h": 4,
                    "w": 6,
                    "x": 12,
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
                    "w": 4,
                    "x": 8,
                    "y": 13,
                    "minH": 3,
                    "minW": 3
                  },
                  "xl": {
                    "h": 4,
                    "w": 6,
                    "x": 12,
                    "y": 5,
                    "minH": 3,
                    "minW": 3
                  },
                  "xxl": {
                    "h": 4,
                    "w": 6,
                    "x": 12,
                    "y": 5,
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
                    "h": 4,
                    "w": 6,
                    "x": 12,
                    "y": 10,
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
                    "w": 4,
                    "x": 8,
                    "y": 16,
                    "minH": 3,
                    "minW": 3
                  },
                  "xl": {
                    "h": 4,
                    "w": 6,
                    "x": 12,
                    "y": 9,
                    "minH": 3,
                    "minW": 3
                  },
                  "xxl": {
                    "h": 4,
                    "w": 6,
                    "x": 12,
                    "y": 9,
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
                    "h": 4,
                    "w": 9,
                    "x": 0,
                    "y": 14,
                    "minH": 3,
                    "minW": 8
                  },
                  "md": {
                    "h": 3,
                    "w": 9,
                    "x": 0,
                    "y": 13,
                    "minH": 3,
                    "minW": 8
                  },
                  "sm": {
                    "h": 3,
                    "w": 8,
                    "x": 0,
                    "y": 13,
                    "minH": 3,
                    "minW": 8
                  },
                  "xl": {
                    "h": 4,
                    "w": 9,
                    "x": 0,
                    "y": 13,
                    "minH": 3,
                    "minW": 8
                  },
                  "xxl": {
                    "h": 3,
                    "w": 9,
                    "x": 0,
                    "y": 13,
                    "minH": 3,
                    "minW": 8
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
                  "lg": {
                    "h": 4,
                    "w": 9,
                    "x": 0,
                    "y": 18,
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
                    "w": 8,
                    "x": 0,
                    "y": 22,
                    "minH": 3,
                    "minW": 8
                  },
                  "xl": {
                    "h": 4,
                    "w": 9,
                    "x": 0,
                    "y": 17,
                    "minH": 3,
                    "minW": 8
                  },
                  "xxl": {
                    "h": 3,
                    "w": 9,
                    "x": 0,
                    "y": 16,
                    "minH": 3,
                    "minW": 8
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
                  "lg": {
                    "h": 4,
                    "w": 9,
                    "x": 9,
                    "y": 14,
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
                    "w": 8,
                    "x": 0,
                    "y": 16,
                    "minH": 3,
                    "minW": 8
                  },
                  "xl": {
                    "h": 4,
                    "w": 9,
                    "x": 9,
                    "y": 13,
                    "minH": 3,
                    "minW": 8
                  },
                  "xxl": {
                    "h": 3,
                    "w": 9,
                    "x": 9,
                    "y": 13,
                    "minH": 3,
                    "minW": 8
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
                  "lg": {
                    "h": 4,
                    "w": 9,
                    "x": 9,
                    "y": 14,
                    "minH": 3,
                    "minW": 8
                  },
                  "md": {
                    "h": 3,
                    "w": 9,
                    "x": 0,
                    "y": 16,
                    "minH": 3,
                    "minW": 8
                  },
                  "sm": {
                    "h": 3,
                    "w": 8,
                    "x": 0,
                    "y": 16,
                    "minH": 3,
                    "minW": 8
                  },
                  "xl": {
                    "h": 4,
                    "w": 9,
                    "x": 9,
                    "y": 17,
                    "minH": 3,
                    "minW": 8
                  },
                  "xxl": {
                    "h": 3,
                    "w": 9,
                    "x": 9,
                    "y": 16,
                    "minH": 3,
                    "minW": 8
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
                  "lg": {
                    "h": 4,
                    "w": 9,
                    "x": 0,
                    "y": 26,
                    "minH": 3,
                    "minW": 8
                  },
                  "md": {
                    "h": 3,
                    "w": 9,
                    "x": 0,
                    "y": 22,
                    "minH": 3,
                    "minW": 8
                  },
                  "sm": {
                    "h": 3,
                    "w": 8,
                    "x": 0,
                    "y": 33,
                    "minH": 3,
                    "minW": 8
                  },
                  "xl": {
                    "h": 4,
                    "w": 9,
                    "x": 0,
                    "y": 21,
                    "minH": 3,
                    "minW": 8
                  },
                  "xxl": {
                    "h": 3,
                    "w": 8,
                    "x": 0,
                    "y": 19,
                    "minH": 3,
                    "minW": 8
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
                  "lg": {
                    "h": 4,
                    "w": 9,
                    "x": 0,
                    "y": 22,
                    "minH": 3,
                    "minW": 8
                  },
                  "md": {
                    "h": 3,
                    "w": 9,
                    "x": 0,
                    "y": 19,
                    "minH": 3,
                    "minW": 8
                  },
                  "sm": {
                    "h": 3,
                    "w": 8,
                    "x": 0,
                    "y": 30,
                    "minH": 3,
                    "minW": 8
                  },
                  "xl": {
                    "h": 4,
                    "w": 9,
                    "x": 0,
                    "y": 25,
                    "minH": 3,
                    "minW": 8
                  },
                  "xxl": {
                    "h": 3,
                    "w": 8,
                    "x": 0,
                    "y": 22,
                    "minH": 3,
                    "minW": 8
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
                  "lg": {
                    "h": 4,
                    "w": 12,
                    "x": 0,
                    "y": 30,
                    "minH": 3,
                    "minW": 8
                  },
                  "md": {
                    "h": 3,
                    "w": 9,
                    "x": 0,
                    "y": 25,
                    "minH": 3,
                    "minW": 8
                  },
                  "sm": {
                    "h": 3,
                    "w": 8,
                    "x": 0,
                    "y": 36,
                    "minH": 3,
                    "minW": 8
                  },
                  "xl": {
                    "h": 4,
                    "w": 9,
                    "x": 9,
                    "y": 29,
                    "minH": 3,
                    "minW": 8
                  },
                  "xxl": {
                    "h": 3,
                    "w": 8,
                    "x": 8,
                    "y": 25,
                    "minH": 3,
                    "minW": 8
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
                  "lg": {
                    "h": 8,
                    "w": 9,
                    "x": 9,
                    "y": 22,
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
                    "w": 8,
                    "x": 9,
                    "y": 25,
                    "minH": 3,
                    "minW": 3
                  },
                  "xl": {
                    "h": 8,
                    "w": 9,
                    "x": 9,
                    "y": 21,
                    "minH": 3,
                    "minW": 3
                  },
                  "xxl": {
                    "h": 6,
                    "w": 10,
                    "x": 8,
                    "y": 19,
                    "minH": 3,
                    "minW": 3
                  }
                },
                "x": {
                  "quantity": "results.properties.catalytic.reactivity.reactants.conversion"
                },
                "y": {
                  "quantity": "results.properties.catalytic.reactivity.products.selectivity"
                },
                "markers": {
                  "color": {
                    "quantity": "results.properties.catalytic.catalyst_characterization.surface_area"
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
                    "y": 34,
                    "minH": 3,
                    "minW": 8
                  },
                  "md": {
                    "h": 3,
                    "w": 9,
                    "x": 0,
                    "y": 28,
                    "minH": 3,
                    "minW": 8
                  },
                  "sm": {
                    "h": 3,
                    "w": 8,
                    "x": 0,
                    "y": 39,
                    "minH": 3,
                    "minW": 8
                  },
                  "xl": {
                    "h": 4,
                    "w": 9,
                    "x": 0,
                    "y": 29,
                    "minH": 3,
                    "minW": 8
                  },
                  "xxl": {
                    "h": 3,
                    "w": 8,
                    "x": 0,
                    "y": 25,
                    "minH": 3,
                    "minW": 8
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
  }
}
