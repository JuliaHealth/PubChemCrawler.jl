# The PUG XML interface is much more painful, but it can run longer queries

"""
    cids = query_substructure_pug(;cid=nothing, smiles=nothing, smarts=nothing,        # specifier for the substructure to search for
                                   maxhits=200_000, poll_interval=10)

Retrieve a list of compounds containing a substructure specified via its `cid`, the SMILES string, or a SMARTS string.

# Example

```julia
julia> using PubChemCrawler

julia> cids = query_substructure_pug(smarts="[r13]Br")   # query brominated 13-atom rings
66-element Vector{Int64}:
  54533707
 153064026
 152829033
...
```

PUG searches can take a while to run (they poll for completion), but conversely they allow more complex, long-running searches
to succeed. See also [`query_substructure`](@ref).
"""
function query_substructure_pug(; kwargs...)
    xdoc = create_substructure_query(; kwargs...)
    cids = submit_substructure_query(xdoc; kwargs...)
    free(xdoc)
    return cids
end

# Implement a common pattern in the XML structure of PUG requests:
# <Tag1>
#   <Tag1_modifierA>
#     <Tag2>
#       <Tag2_modifierB>
# ...
function nest2(parent, tag, mod)
    nodetag = new_child(parent, tag)
    child = new_child(nodetag, tag*'_'*mod)
    return nodetag, child
end
function nest2(parent, tag, mods...)
    nodetag = new_child(parent, tag)
    children = [new_child(nodetag, tag*'_'*mod) for mod in mods]
    return nodetag, children
end



function create_substructure_query(;cid=nothing, smiles=nothing, smarts=nothing,          # inputs
    stereo="ignore",isotopes::Bool=false,charges::Bool=false,tautomers::Bool=false,rings::Bool=false,bonds::Bool=true,chains::Bool=true,hydrogen::Bool=false,
    maxhits::Int=2_000_000,
    kwargs...)

    query_data = if cid !== nothing
        (smiles === nothing && smarts === nothing) || throw(ArgumentError("only one of cid, smiles, or smarts can be specified"))
        string(cid)
    elseif smiles !== nothing
        smarts === nothing || throw(ArgumentError("only one of cid, smiles, or smarts can be specified"))
        smiles
    elseif smarts !== nothing
        smarts
    else
        throw(ArgumentError("one of cid, smiles, or smarts must be specified"))
    end
    stereoidx = 0
    if stereo != "ignore"
        stereoidx = findfirst(isequal(stereo), ["exact", "relative", "non-conflicting"])
        stereoidx === nothing && throw(ArgumentError("""stereo must be one of "ignore", "exact", "relative", or "non-conflicting", got $stereo"""))
    end

    xdoc = XMLDocument()
    xroot = create_root(xdoc, "PCT-Data")
    rootinput = new_child(xroot, "PCT-Data_input")
    input, inputquery = nest2(rootinput, "PCT-InputData", "query")
    query, querytype = nest2(inputquery, "PCT-Query", "type")
    querytype2, querytype2css = nest2(querytype, "PCT-QueryType", "css")
    querycomp, fields = nest2(querytype2css, "PCT-QueryCompoundCS", "query", "type", "results")
    add_text(new_child(fields[1], "PCT-QueryCompoundCS_query_data"), query_data)  # under "query"
    subss = new_child(fields[2], "PCT-QueryCompoundCS_type_subss")                # under "type"

    ssstr = "PCT-CSStructure"
    substruct = new_child(subss, ssstr)
    if stereoidx != 0
        cstereo = new_child(substruct, ssstr*"_stereo")
        set_attribute(cstereo, "value", stereo)
        add_text(cstereo, string(stereoidx))
    end
    for (name, val) in (("isotopes", isotopes), ("charges", charges), ("tautomers", tautomers), ("rings", rings), ("bonds", bonds), ("chains", chains), ("hydrogen", hydrogen))
        child = new_child(substruct, ssstr*'_'*name)
        set_attribute(child, "value", repr(val))
    end

    add_text(fields[3], string(maxhits))               # under "results"

    return xdoc
end

function create_poll(id)
    xdoc = XMLDocument()
    xroot = create_root(xdoc, "PCT-Data")
    rootinput = new_child(xroot, "PCT-Data_input")
    input, inputreq = nest2(rootinput, "PCT-InputData", "request")
    req, reqdata = nest2(inputreq, "PCT-Request", "reqid", "type")
    add_text(reqdata[1], id)
    set_attribute(reqdata[2], "value", "status")
    return xdoc
end

# function create_download(entrez; format="CSV", compression="gzip")
#     xdoc = XMLDocument()
#     xroot = create_root(xdoc, "PCT-Data")
#     rootinput = new_child(xroot, "PCT-Data_input")
#     input, inputreq = nest2(rootinput, "PCT-InputData", "download")
#     req, reqdata = nest2(inputreq, "PCT-Download", "uids", "format", "compression")
#     query, centrez = nest2(reqdata[1], "PCT-QueryUids", "entrez")
#     _, fields = nest2(centrez, "PCT-Entrez", "db", "query-key", "webenv")
#     add_text(fields[1], content(entrez["PCT-Entrez_db"][1]))
#     add_text(fields[2], content(entrez["PCT-Entrez_query-key"][1]))
#     add_text(fields[3], content(entrez["PCT-Entrez_webenv"][1]))
#     set_attribute(reqdata[2], "value", format)
#     set_attribute(reqdata[3], "value", compression)
#     return xdoc
# end

function get_top(xdoc)
    xdoc_root = root(xdoc)
    return xdoc_root["PCT-Data_output"][1]["PCT-OutputData"][1]
end
function get_status(top)
    elstatus = top["PCT-OutputData_status"][1]["PCT-Status-Message"][1]["PCT-Status-Message_status"][1]["PCT-Status"][1]
    return attribute(elstatus, "value")
end
function get_id(top)
    return content(top["PCT-OutputData_output"][1]["PCT-OutputData_output_waiting"][1]["PCT-Waiting"][1]["PCT-Waiting_reqid"][1])
end
function get_entrez(top)
    return top["PCT-OutputData_output"][1]["PCT-OutputData_output_entrez"][1]["PCT-Entrez"][1]
end

function submit_substructure_query(xdoc; poll_interval=10, kwargs...)
    r = HTTP.request("POST", "https://pubchem.ncbi.nlm.nih.gov/pug/pug.cgi", [], string(xdoc); kwargs...)
    xresp = parse_string(String(r.body))
    top = get_top(xresp)
    status = get_status(top)
    if status == "queued"
        # Wait for the results
        id = get_id(top)
        pollreq = create_poll(id)
        while true
            free(xresp)
            sleep(poll_interval)
            # poll for completion
            r = HTTP.request("POST", "https://pubchem.ncbi.nlm.nih.gov/pug/pug.cgi", [], string(pollreq); kwargs...)
            xresp = parse_string(String(r.body))
            top = get_top(xresp)
            status = get_status(top)
            status != "queued" && status != "running" && break
        end
        free(pollreq)
    end
    if status == "hit-limit"
        @warn "maxhits was hit, results are partial"
    elseif status != "success"
        error(string(xresp))
    end
    entrez = get_entrez(top)
    # Now retrieve the results
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?retmode=text&rettype=uilist&WebEnvRq=1&db=pccompound"
    url *= "&query_key=" * content(entrez["PCT-Entrez_query-key"][1])
    url *= "&WebEnv=" * content(entrez["PCT-Entrez_webenv"][1])
    free(xresp)
    r = HTTP.request("GET", url)
    return parse.(Int, split(chomp(String(r.body)), '\n'))
    # xdl = create_download(entrez; format=format)
    # r = HTTP.request("POST", "https://pubchem.ncbi.nlm.nih.gov/pug/pug.cgi", [], string(xdl); kwargs...)
    # result = String(r.body)
    # free(xdl)
end
