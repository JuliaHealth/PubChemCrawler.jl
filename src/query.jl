const prolog = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/"

"""
    get_cids(; name=nothing, smiles=nothing, cas_number=nothing,kwargs...)

Return all the PubChem **c**ompound **id**entification numbers for the specified compound.

- `get_cid` returns a single identifier and fails if there are multiple results.
- `get_cids` returns a vector of identifiers, containing all the identifiers that match

Queries on `cas_number` often return multiple `cids`.

Examples:
```
julia> get_cids(name="2-nonenal")
3-element Vector{Int64}:
 5283335
   17166
 5354833

julia> get_cid(name="2-nonenal")
ERROR: ArgumentError: Collection has multiple elements, must contain exactly 1 element

julia> get_cids(cas_number="50-78-2")
4-element Vector{Int64}:
     2244
    67252
  3434975
 12280114

```
"""
function get_cids(; name=nothing, smiles=nothing, cas_number=nothing,
                   kwargs...)
    input = "compound/"
    name !== nothing && (input *= "name/$(HTTP.escapeuri(name))/")
    smiles !== nothing && (input *= "smiles/$((smiles))/")
    cas_number !== nothing && (input *= "xref/RN/$(cas_number)/")
    url = prolog * input * "cids/TXT"
    r = HTTP.request("GET", url; kwargs...)
    return parse.(Int,split(chomp(String(r.body)), '\n'))
end

"""
    get_cid(; name=nothing, smiles=nothing, cas_number=nothing, kwargs...)

Return the PubChem **c**ompound **id**entification number for the specified compound.

Examples:
```
julia> cid = get_cid(name="glucose")
5793

julia> cid = get_cid(smiles="C([C@@H]1[C@H]([C@@H]([C@H](C(O1)O)O)O)O)O")
5793
```
"""
get_cid = only ∘ get_cids


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

# Example

```
julia> using PubChemCrawler, CSV, DataFrames

julia> cid = get_cid(name="estriol")
5756

julia> df = CSV.File(query_substructure(;cid)) |> DataFrame      # on Julia 1.0, use `(;cid=cid)`
11607×4 DataFrame
│ Row  │ CID       │ MolecularFormula │ MolecularWeight │ XLogP    │
│      │ $Int     │ String           │ Float64         │ Float64? │
├──────┼───────────┼──────────────────┼─────────────────┼──────────┤
│ 1    │ 5756      │ C18H24O3         │ 288.4           │ 2.5      │
│ 2    │ 5281904   │ C24H32O9         │ 464.5           │ 1.1      │
│ 3    │ 27125     │ C18H24O4         │ 304.4           │ 1.5      │
...
```

will query for derivatives of [estriol](https://en.wikipedia.org/wiki/Estriol).

!!! info
    For complex queries that risk timing out, consider [`query_substructure_pug`](@ref) in combination with [`get_for_cids`](@ref).
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
    msg = get_for_cids(cids; properties|xrefs|cids_type|record_type, output="CSV")

Retrieve the given `properties`, `xrefs`, CIDs, or records, respectively, for a list of compounds specified by their `cids`.
The documentation for these traits can be found at http://pubchemdocs.ncbi.nlm.nih.gov/pug-rest; this URL will be referred to as
PUGREST below.

- `properties` include structural features like the molecular formula, number of undefined stereocenters, and so on.
  Specify these as a comma-separated list from among the choices in PUGREST under "Compound Property Tables".
- `xrefs` ("cross-references") include identifiers used by other databases, e.g., the CAS (Registry) number, PubMedID, and so on.
  The supported values for `xrefs` are available at PUGREST under "XRefs".
- `cids_type` is used to retrieve CIDs for compounds related to those specified in `cids`; see PUGREST under "SIDS / CIDS / AIDS".
- `record_type` is used to retrieve data files and to specify options for these files, e.g., 2d or 3d SDF files.
  See PUGREST under "Full-record Retrieval".

`output` specifies the output format. Not all options are applicable to all queries; for example, "CSV" is appropriate for
`properties` queries but "SDF" might be used for a `record_type` query. See PUGREST, "Output".

# Examples

```
julia> using PubChemCrawler, CSV, DataFrames, JSON3

julia> cids = [get_cid(name="cyclic guanosine monophosphate"), get_cid(name="aspirin")]
2-element Array{$Int,1}:
 135398570
      2244

julia> CSV.File(get_for_cids(cids; properties="MolecularFormula,XLogP", output="CSV")) |> DataFrame
2×3 DataFrame
 Row │ CID        MolecularFormula  XLogP
     │ $Int      String            Float64
─────┼──────────────────────────────────────
   1 │ 135398570  C10H12N5O7P          -3.4
   2 │      2244  C9H8O4                1.2

julia> open("/tmp/aspirin_3d.sdf", "w") do io    # save the 3d SDF file for aspirin (CID 2244)
           write(io, get_for_cids(2244; record_type="3d", output="SDF"))
       end
4055

julia> dct = JSON3.read(get_for_cids(cids; xrefs="RN,", output="JSON"));   # get the Registry Number(s) (CAS)

julia> dct[:InformationList][:Information]
2-element JSON3.Array{JSON3.Object,Array{UInt8,1},SubArray{$UInt,1,Array{$UInt,1},Tuple{UnitRange{$Int}},true}}:
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
                      cids_type=nothing,
                      record_type=nothing,
                      output="CSV",
                      kwargs...)
    url = prolog * "compound/cid/"
    if xrefs === nothing
        if properties !== nothing
            url *= canonicalize_properties("property/" * properties)
        elseif cids_type !== nothing
            url *= "cids/"
        end
    else
        properties === nothing || error("cannot specify both xref and properties in a single query")
        url *= canonicalize_properties("xrefs/" * xrefs)
    end
    url = joinpath(url, output)
    if cids_type !== nothing
        url *= "?cids_type=" * cids_type
    end
    if record_type !== nothing
        url *= "?record_type=" * record_type
    end
    r = HTTP.request("POST", url, ["Content-Type"=>"application/x-www-form-urlencoded"], "cid="*join(string.(cids), ","); kwargs...)
    return r.body
end

get_for_cids(cid::Int; kwargs...) = get_for_cids([cid]; kwargs...)

"""
   `pug(args...; silent = true, escape_args = true, return_text = true, status_exception = false, kwargs...)`

Generate a PUG endpoint and call it. The details about PUG endpoints are described here: <https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest>

Keyword arguments:

- `escape_args = true`, URL encodes each argument before generating the endpoint.
Setting this false is useful when copy-pasting an existing PUG endpoint, e.g. from documentation.

- `silent = false` print the pug URL called.

-  `return_text = true`, call `String` on the output to return a string rather than a byte vector.

-  `status_exception = false`, tell HTTP.jl to not throw an exception on return codes >= 300.

Other keyword arguments are passed on to `HTTP.request`.

Examples:

```
julia> pug(:compound, :name, "ethanol", :cids, :txt, silent = false, return_text = true)
[ Info: https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/ethanol/cids/txt
"702"

julia> pug("compound/cid/2244", :cids, :txt, escape_args = false, silent = false, return_text = true)
[ Info: https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/2244/cids/txt
"2244"

julia> pug(:compound, :smiles, "C([C@@H]1[C@H]([C@@H]([C@H](C(O1)O)O)O)O)O", :cids, :txt, return_text = true)
"5793"

julia> pug(:compound, :cid, 708, :txt, return_text = true, status_exception = false)
"Status: 400\nCode: PUGREST.BadRequest\nMessage: Invalid output format\nDetail: Full-record output format must be one of ASNT/B, XML, JSON(P), SDF, or PNG"
```
"""
function pug(args...; silent = true, escape_args = true, return_text = false, status_exception = true, kwargs...)
    args =  replace.(string.(args), r"/$" => "", r"^/" => "")
    escape_args && (args = HTTP.escapeuri.(args))
    pug_string = join(string.(args), "/")
    url = prolog * pug_string
    silent || @info url
    r = HTTP.request("GET", url; status_exception, kwargs...)
    b = return_text ? chomp(String(r.body)) : r.body
    return b
end

"""
    synonyms = get_synonyms(name="glucose")
    synonyms = get_synonyms(smiles="C([C@@H]1[C@H]([C@@H]([C@H](C(O1)O)O)O)O)O")
    synonyms = get_synonyms(cid=5793)

Return a list of substance or compound synonyms.
"""
function get_synonyms(; name=nothing, cid=nothing, smiles=nothing,
         kwargs...)::Vector{String}
    # inputs
    input = "compound/"
    if name !== nothing
        (smiles === nothing && cid === nothing) || throw(ArgumentError("only one of name, cid, or smiles can be specified"))
        input *= "name/$(HTTP.escapeuri(name))/"
    elseif cid !== nothing
        (smiles === nothing) || throw(ArgumentError("only one of name, cid, or smiles can be specified"))
        input *= "cid/$cid/"
    elseif smiles !== nothing
        input *= "smiles/$smiles/"
    else
        throw(ArgumentError("one of name, cid, or smiles must be specified"))
    end
    url = prolog * input * "synonyms/TXT"
    r = HTTP.request("GET", url; kwargs...)
    return split(chomp(String(r.body)), "\n")
end
