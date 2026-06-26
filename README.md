# PubChemCrawler.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaHealth.github.io/PubChemCrawler.jl/stable)
[![Build Status](https://github.com/JuliaHealth/PubChemCrawler.jl/workflows/CI/badge.svg)](https://github.com/JuliaHealth/PubChemCrawler.jl/actions)
[![Coverage](https://codecov.io/gh/JuliaHealth/PubChemCrawler.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaHealth/PubChemCrawler.jl)

PubChemCrawler.jl makes it easy to query the [PubChem database](https://pubchem.ncbi.nlm.nih.gov/) from Julia. This package provides a simple interface to search for chemical compounds, retrieve their properties, and perform substructure searches.

## Features

- **Compound Lookup**: Search for compounds by name, SMILES, or CAS number
- **Property Retrieval**: Get molecular formulas, weights, and other properties
- **Substructure Search**: Find compounds containing specific molecular substructures
- **Synonym Retrieval**: Access chemical names and synonyms
- **Multiple Output Formats**: CSV, JSON, TXT, SDF, and more

## Installation

Install PubChemCrawler.jl using Julia's package manager:

```julia
using Pkg
Pkg.add("PubChemCrawler")
```

Or in the Julia REPL, press `]` to enter package mode and run:

```
add PubChemCrawler
```

## Quick Start

```julia
using PubChemCrawler

# Get compound ID for a chemical by name
cid = get_cid(name="aspirin")

# Get compound ID using SMILES notation
cid = get_cid(smiles="CC(=O)OC1=CC=CC=C1C(=O)O")

# Get all compound IDs for a CAS number
cids = get_cids(cas_number="50-78-2")
```

## Usage Examples

For detailed usage examples including:
- Finding compound IDs
- Retrieving compound properties
- Substructure searching
- Getting synonyms
- Parsing chemical formulas

See the [Getting Started guide](https://juliahealth.org/PubChemCrawler.jl/stable/#Getting-started).

## Documentation

For detailed documentation, including all available functions and options, visit the [official documentation](https://JuliaHealth.github.io/PubChemCrawler.jl/stable).

## Contributing

We welcome contributions! Please see our [Contributing Guidelines](CONTRIBUTING.md) for details on how to:

- Report issues
- Submit pull requests
- Follow code style guidelines (Blue style)

## License

This package is licensed under the MIT License. See [LICENSE](LICENSE) for details.

