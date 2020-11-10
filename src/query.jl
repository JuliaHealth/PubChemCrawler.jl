const prolog = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/"

function get_cid(; name=nothing)                   # inputs
    input = "compound/"
    name !== nothing && (input *= "name/$name/")
    url = prolog * input * "cids/TXT"
    r = HTTP.request("GET", url)
    return parse(Int, chomp(String(r.body)))
end

function query_substructure(;cid=nothing, smiles=nothing, name=nothing,   # inputs
                             properties="MolecularFormula,MolecularWeight,XLogP,", # http://pubchemdocs.ncbi.nlm.nih.gov/pug-rest, "Compound Property Tables"
                             output="CSV")
    input = "compound/fastsubstructure/"
    if cid !== nothing
        input *= "cid/$cid/"
    elseif smiles !== nothing
        input *= "smiles/$smiles/"
    else
        error("must specify input method, cid, smiles, or name")
    end
    props = canonicalize_properties("property/" * properties)
    url = prolog * input * props * output * "?StripHydrogen=true"
    r = HTTP.request("GET", url)
    return r.body
end

function get_for_cids(cids;
                      properties="CanonicalSMILES,IsomericSMILES,",
                      output="CSV")
    input = "compound/cid/"
    props = canonicalize_properties("property/" * properties)
    url = prolog * input * props * output
    r = HTTP.request("POST", url, ["Content-Type"=>"application/x-www-form-urlencoded"], "cid="*join(string.(cids), ","))
    return r.body
end
