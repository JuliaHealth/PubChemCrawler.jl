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
- Requests made via REST time out after 30s. The PUG XML interface does not have this limitation. For substructure searches, [`query_substructure_pug`](@ref) is recommended.

## Getting started

Some queries make use of the CID, the **C**ompound **ID**entifier, which you can obtain in a variety of ways.
Let's get the CID for aspirin:

```julia
julia> cid = get_cid(name="aspirin")
2244

julia> cid = get_cid(smiles="CC(=O)OC1=CC=CC=C1C(=O)O")   # use the SMILES string
2244
```

You can then retrieve individual properties:

```julia
julia> smiles = chomp(String(get_for_cids(2244, properties="CanonicalSMILES", output="TXT")))
"CC(=O)OC1=CC=CC=C1C(=O)O"
```

or a list of properties:

```julia
julia> using CSV, DataFrames

julia> df = CSV.File(get_for_cids(2244; properties="MolecularFormula,MolecularWeight,XLogP,IsomericSMILES", output="CSV")) |> DataFrame
1×5 DataFrame
│ Row │ CID   │ MolecularFormula │ MolecularWeight │ XLogP   │ IsomericSMILES           │
│     │ $Int │ String           │ Float64         │ Float64 │ String                   │
├─────┼───────┼──────────────────┼─────────────────┼─────────┼──────────────────────────┤
│ 1   │ 2244  │ C9H8O4           │ 180.16          │ 1.2     │ CC(=O)OC1=CC=CC=C1C(=O)O │
```

You can query properties for a whole list of `cids`.

If your query returns multiple `cids`, you need to use `get_cids`:

``` julia
julia> cids = get_cids(cas_number="50-78-2")
4-element Vector{Int64}:
     2244
    67252
  3434975
 12280114
```

You can also download structure data and save it to a file. This saves a 3d conformer for aspirin:

```julia
julia> open("/tmp/aspirin.sdf", "w") do io
           write(io, get_for_cids(2244, output="SDF", record_type="3d"))
       end
3637
```

Finally, you can perform substructure searches. Let's retrieve up to 10 [bicyclic](https://en.wikipedia.org/wiki/Bicyclic_molecule) compounds using a [SMARTS](https://en.wikipedia.org/wiki/SMILES_arbitrary_target_specification) search:

```julia
julia> cids = query_substructure_pug(smarts = "[\$([*R2]([*R])([*R])([*R]))].[\$([*R2]([*R])([*R])([*R]))]", maxhits = 10)
┌ Warning: maxhits was hit, results are partial
└ @ PubChemCrawler ~/.julia/dev/PubChemCrawler/src/pugxml.jl:164
10-element Vector{$Int}:
 135398658
   5280795
      5430
      5143
  54675779
   5280961
   5280804
   5280793
   5280343
   3034034
```

Note that Julia (not this package) requires the SMARTS string characters `$` be escaped.

## API

### Queries

```@docs
get_cid
get_cids
query_substructure_pug
query_substructure
get_for_cids
pug
get_synonyms
```

### Utilities

```@docs
parse_formula
atomregex
```
