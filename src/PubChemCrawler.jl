module PubChemCrawler

using HTTP

export atomregex, parse_formula, get_cid, get_for_cids, query_substructure

include("utils.jl")
include("query.jl")

end
