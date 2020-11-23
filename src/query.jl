const prolog = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/"

"""
    cid = get_cid(name="glucose")
    cid = get_cid(smiles="C([C@@H]1[C@H]([C@@H]([C@H](C(O1)O)O)O)O)O")

Return the PubChem **c**ompound **id**entification number for the specified compound.
"""
function get_cid(; name=nothing, smiles=nothing,                   # inputs
                   kwargs...)
    input = "compound/"
    name !== nothing && (input *= "name/$(HTTP.escapeuri(name))/")
    smiles !== nothing && (input *= "smiles/$((smiles))/")
    url = prolog * input * "cids/TXT"
    r = HTTP.request("GET", url; kwargs...)
    return parse(Int, chomp(String(r.body)))
end

"""
    msg = query_substructure(;cid=nothing, smiles=nothing, smarts=nothing,           # specifier for the substructure to search for
                              properties="MolecularFormula,MolecularWeight,XLogP,",  # properties to retrieve
                              output="CSV")                                          # output format

Perform a substructure search of the entire database. You can specify the target via its `cid`, the SMILES string, or a SMARTS string.
Specify the `properties` you want to retrieve as a comma-separated list from among the choices in
http://pubchemdocs.ncbi.nlm.nih.gov/pug-rest, "Compound Property Tables".
Requesting more properties takes more time.

The output is a `Vector{UInt8}`. For `output="CSV"`, a good choice to generate a manipulable result is
`DataFrame(CSV.File(msg))` from the DataFrames and CSV packages, respectively.
Alternatively `String(msg)` will convert it to a string, which you can write to a file.

For complex queries that risk timing out, consider [`query_substructure_pug`](@ref).

# Example

```
julia> using PubChemCrawler, CSV, DataFrames

julia> cid = get_cid(name="estriol")
5756

julia> df = CSV.File(query_substructure(;cid)) |> DataFrame      # on Julia 1.0, use `(;cid=cid)`
11607×4 DataFrame
│ Row  │ CID       │ MolecularFormula │ MolecularWeight │ XLogP    │
│      │ Int64     │ String           │ Float64         │ Float64? │
├──────┼───────────┼──────────────────┼─────────────────┼──────────┤
│ 1    │ 5756      │ C18H24O3         │ 288.4           │ 2.5      │
│ 2    │ 5281904   │ C24H32O9         │ 464.5           │ 1.1      │
│ 3    │ 27125     │ C18H24O4         │ 304.4           │ 1.5      │
...
```

will query for derivatives of [estriol](https://en.wikipedia.org/wiki/Estriol).
"""
function query_substructure(;cid=nothing, smiles=nothing, smarts=nothing,          # inputs
                             properties="MolecularFormula,MolecularWeight,XLogP,", # http://pubchemdocs.ncbi.nlm.nih.gov/pug-rest, "Compound Property Tables"
                             output="CSV",
                             kwargs...)
    input = "compound/fastsubstructure/"
    if cid !== nothing
        (smiles === nothing && smarts === nothing) || throw(ArgumentError("only one of cid, smiles, or smarts can be specified"))
        input *= "cid/$cid/"
    elseif smiles !== nothing
        smarts === nothing || throw(ArgumentError("only one of cid, smiles, or smarts can be specified"))
        input *= "smiles/$smiles/"
    elseif smarts !== nothing
        input *= "smarts/$smarts/"
    else
        throw(ArgumentError("one of cid, smiles, or smarts must be specified"))
    end
    props = canonicalize_properties("property/" * properties)
    url = prolog * input * props * output * "?StripHydrogen=true"
    r = HTTP.request("GET", url; kwargs...)
    return r.body
end

"""
    msg = get_for_cids(cids; properties=nothing, xrefs=nothing, output="CSV")

Retrieve the given `properties` or `xrefs` for a list of compounds specified by their `cids`.

See [`query_substructure`](@ref) for information about the arguments and return value.
The supported values for `xrefs` are available at https://pubchemdocs.ncbi.nlm.nih.gov/pug-rest under "XRefs".

# Example

```
julia> using PubChem, JSON3

julia> cids = [get_cid(name="cyclic guanosine monophosphate"), get_cid(name="aspirin")]
2-element Array{Int64,1}:
 135398570
      2244

julia> dct = JSON3.read(get_for_cids(cids; xrefs="RN,", output="JSON"))   # get the Registry Number(s) (CAS)
JSON3.Object{Array{UInt8,1},Array{UInt64,1}} with 1 entry:
  :InformationList => {…

julia> dct[:InformationList][:Information]
2-element JSON3.Array{JSON3.Object,Array{UInt8,1},SubArray{UInt64,1,Array{UInt64,1},Tuple{UnitRange{Int64}},true}}:
 {
   "CID": 135398570,
    "RN": [
            "40732-48-7",
            "7665-99-8"
          ]
}
 {
   "CID": 2244,
    "RN": [
            "11126-35-5",
            "156865-15-5",
            "50-78-2",
            "52080-78-1",
            "921943-73-9",
            "98201-60-6",
            "99512-66-0"
          ]
}
```
"""
function get_for_cids(cids;
                      properties=nothing,
                      xrefs=nothing,
                      output="CSV",
                      record_type=nothing,
                      kwargs...)
    url = prolog * "compound/cid/"
    if xrefs === nothing
        if properties !== nothing
            url *= canonicalize_properties("property/" * properties)
        end
    else
        properties === nothing || error("cannot specify both xref and properties in a single query")
        url *= canonicalize_properties("xrefs/" * xrefs)
    end
    url = joinpath(url, output)
    if record_type !== nothing
        url *= "?record_type=" * record_type
    end
    r = HTTP.request("POST", url, ["Content-Type"=>"application/x-www-form-urlencoded"], "cid="*join(string.(cids), ","); kwargs...)
    return r.body
end

get_for_cids(cid::Int; kwargs...) = get_for_cids([cid]; kwargs...)
