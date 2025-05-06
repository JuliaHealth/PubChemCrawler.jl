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
    sdfs = get_conformers_for_cid(cid, n = Inf)
Retrieve 3D records for up to `n` conformers for a compound specified by its `cid`. Conformer ordering is such that the first
"n" conformers selected represent the overall diversity of the conformer model for a compound. A description of PubChem's diverse 
conformer ordering can be found at https://pubchem.ncbi.nlm.nih.gov/release3d.html. 
# Example
````
julia> sdfs = get_conformers_for_cid(get_cid(name="aspirin"), 5);   # get data for the first 5 conformers for aspirin
julia> for (i,sdf) in enumerate(sdfs)   # save the 3d SDF files for each retrieved conformer
           open("tmp/aspirin_conf"*string(i)*"_3d.sdf", "w") do io
               write(io, sdf)
           end
       end
````
"""
function get_conformers_for_cid(cid, n = Inf)
    url = prolog * "compound/cid/" * string(cid) * "/conformers/XML"
    r = HTTP.request("GET", url)
    xdoc = parse_string(String(r.body)) 
    xroot = root(xdoc)
    confids = []
    for c in child_nodes(xroot)
        if is_elementnode(c)
            e = XMLElement(c) 
            els = get_elements_by_tagname(e, "ConformerID")
            for (i,el) in enumerate(els)
                i > n && break
                push!(confids,content(el))
            end
        end
    end
    confsdfs = []
    for confid in confids
        confurl = prolog * "conformers/" * string(confid) * "/SDF"
        confr = HTTP.request("GET", confurl)
        push!(confsdfs, confr.body)
    end
    return confsdfs
end
