```@meta
CurrentModule = PubChemCrawler
```

# PubChemCrawler

PubChemCrawler makes it easier to search the [PubChem database](https://pubchem.ncbi.nlm.nih.gov/) from Julia.
You can use it to access information about particular compounds or query substructures.

The package supports only a subset of the available [functionality](https://pubchemdocs.ncbi.nlm.nih.gov/pug-rest),
but it is fairly straightforward to expand to other types of query. Pull requests are welcome!
If you do want to make improvements to this package, this [tutorial](https://pubchemdocs.ncbi.nlm.nih.gov/pug-rest-tutorial)
might help you get started.

## Before you start: be aware of search limits

PubChem places significant limits on requests:

- No more than 5 requests per second
- No more than 400 requests per minute
- No longer than 300 second running time per minute
- Requests made via REST (the API targeted by this package) time out after 30s. There seems to be some caching, though, so sometimes it works to just resubmit the request.

The latter constraint is worth a bit of thought. Substructure search requests can sometimes generate hundreds of thousands of hits,
and asking for lots of different properties for many different compounds can prolong the request beyond the timeout limit.
One recommendation is to pursue a "staged" approach, where you first request just a few simple properties
(e.g., `MolecularFormula` rather than the full `SMILES` string) and use that information to narrow the list of candidates.
With a smaller list, it becomes feasible to request more detailed information.

## Getting started

Some queries make use of the CID, the **C**ompound **ID**entifier, which you can obtain in a variety of ways.
Let's get the CID for aspirin:

```julia
julia> cid = get_cid(name="aspirin")
2244

julia> cid = get_cid(smiles="CC(=O)OC1=CC=CC=C1C(=O)O")   # use the SMILES string
2244
```

You can then retrieve a list of other properties:

```julia
julia> df = CSV.File(get_for_cids(cid; properties="MolecularFormula,MolecularWeight,XLogP,IsomericSMILES", output="CSV")) |> DataFrame
1×5 DataFrame
│ Row │ CID   │ MolecularFormula │ MolecularWeight │ XLogP   │ IsomericSMILES           │
│     │ Int64 │ String           │ Float64         │ Float64 │ String                   │
├─────┼───────┼──────────────────┼─────────────────┼─────────┼──────────────────────────┤
│ 1   │ 2244  │ C9H8O4           │ 180.16          │ 1.2     │ CC(=O)OC1=CC=CC=C1C(=O)O │
```

You can query such properties for a list of `cids`.

## API

### Queries

```@docs
get_cid
query_substructure
get_for_cids
```

### Utilities

```@docs
parse_formula
atomregex
```
