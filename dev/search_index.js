var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = PubChemCrawler","category":"page"},{"location":"#PubChemCrawler","page":"Home","title":"PubChemCrawler","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"PubChemCrawler makes it easier to search the PubChem database from Julia. You can use it to access information about particular compounds or query substructures.","category":"page"},{"location":"","page":"Home","title":"Home","text":"The package supports only a subset of the available functionality, but it is fairly straightforward to expand to other types of query. Pull requests are welcome! If you do want to make improvements to this package, this tutorial might help you get started.","category":"page"},{"location":"#Before-you-start:-be-aware-of-search-limits","page":"Home","title":"Before you start: be aware of search limits","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"PubChem places significant limits on requests:","category":"page"},{"location":"","page":"Home","title":"Home","text":"No more than 5 requests per second\nNo more than 400 requests per minute\nNo longer than 300 second running time per minute\nRequests made via REST time out after 30s. The PUG XML interface does not have this limitation. For substructure searches, query_substructure_pug is recommended.","category":"page"},{"location":"#Getting-started","page":"Home","title":"Getting started","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Some queries make use of the CID, the Compound IDentifier, which you can obtain in a variety of ways. Let's get the CID for aspirin:","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> cid = get_cid(name=\"aspirin\")\n2244\n\njulia> cid = get_cid(smiles=\"CC(=O)OC1=CC=CC=C1C(=O)O\")   # use the SMILES string\n2244","category":"page"},{"location":"","page":"Home","title":"Home","text":"You can then retrieve individual properties:","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> smiles = chomp(String(get_for_cids(2244, properties=\"CanonicalSMILES\", output=\"TXT\")))\n\"CC(=O)OC1=CC=CC=C1C(=O)O\"","category":"page"},{"location":"","page":"Home","title":"Home","text":"or a list of properties:","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> using CSV, DataFrames\n\njulia> df = CSV.File(get_for_cids(2244; properties=\"MolecularFormula,MolecularWeight,XLogP,IsomericSMILES\", output=\"CSV\")) |> DataFrame\n1×5 DataFrame\n│ Row │ CID   │ MolecularFormula │ MolecularWeight │ XLogP   │ IsomericSMILES           │\n│     │ $Int │ String           │ Float64         │ Float64 │ String                   │\n├─────┼───────┼──────────────────┼─────────────────┼─────────┼──────────────────────────┤\n│ 1   │ 2244  │ C9H8O4           │ 180.16          │ 1.2     │ CC(=O)OC1=CC=CC=C1C(=O)O │","category":"page"},{"location":"","page":"Home","title":"Home","text":"You can query properties for a whole list of cids.","category":"page"},{"location":"","page":"Home","title":"Home","text":"You can also download structure data and save it to a file. This saves a 3d conformer for aspirin:","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> open(\"/tmp/aspirin.sdf\", \"w\") do io\n           write(io, get_for_cids(2244, output=\"SDF\", record_type=\"3d\"))\n       end\n3637","category":"page"},{"location":"","page":"Home","title":"Home","text":"Finally, you can perform substructure searches. Let's retrieve up to 10 bicyclic compounds using a SMARTS search:","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> cids = query_substructure_pug(smarts = \"[\\$([*R2]([*R])([*R])([*R]))].[\\$([*R2]([*R])([*R])([*R]))]\", maxhits = 10)\n┌ Warning: maxhits was hit, results are partial\n└ @ PubChemCrawler ~/.julia/dev/PubChemCrawler/src/pugxml.jl:164\n10-element Vector{$Int}:\n 135398658\n   5280795\n      5430\n      5143\n  54675779\n   5280961\n   5280804\n   5280793\n   5280343\n   3034034","category":"page"},{"location":"","page":"Home","title":"Home","text":"Note that Julia (not this package) requires the SMARTS string characters $ be escaped.","category":"page"},{"location":"#API","page":"Home","title":"API","text":"","category":"section"},{"location":"#Queries","page":"Home","title":"Queries","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"get_cid\nquery_substructure_pug\nquery_substructure\nget_for_cids","category":"page"},{"location":"#PubChemCrawler.get_cid","page":"Home","title":"PubChemCrawler.get_cid","text":"cid = get_cid(name=\"glucose\")\ncid = get_cid(smiles=\"C([C@@H]1[C@H]([C@@H]([C@H](C(O1)O)O)O)O)O\")\n\nReturn the PubChem compound identification number for the specified compound.\n\n\n\n\n\n","category":"function"},{"location":"#PubChemCrawler.query_substructure_pug","page":"Home","title":"PubChemCrawler.query_substructure_pug","text":"cids = query_substructure_pug(;cid=nothing, smiles=nothing, smarts=nothing,        # specifier for the substructure to search for\n                               maxhits=200_000, poll_interval=10)\n\nRetrieve a list of compounds containing a substructure specified via its cid, the SMILES string, or a SMARTS string.\n\nExample\n\njulia> using PubChemCrawler\n\njulia> cids = query_substructure_pug(smarts=\"[r13]Br\")   # query brominated 13-atom rings\n66-element Vector{Int64}:\n  54533707\n 153064026\n 152829033\n...\n\nPUG searches can take a while to run (they poll for completion), but conversely they allow more complex, long-running searches to succeed. See also query_substructure.\n\n\n\n\n\n","category":"function"},{"location":"#PubChemCrawler.query_substructure","page":"Home","title":"PubChemCrawler.query_substructure","text":"msg = query_substructure(;cid=nothing, smiles=nothing, smarts=nothing,           # specifier for the substructure to search for\n                          properties=\"MolecularFormula,MolecularWeight,XLogP,\",  # properties to retrieve\n                          output=\"CSV\")                                          # output format\n\nPerform a substructure search of the entire database. You can specify the target via its cid, the SMILES string, or a SMARTS string. Specify the properties you want to retrieve as a comma-separated list from among the choices in http://pubchemdocs.ncbi.nlm.nih.gov/pug-rest, \"Compound Property Tables\". Requesting more properties takes more time.\n\nThe output is a Vector{UInt8}. For output=\"CSV\", a good choice to generate a manipulable result is DataFrame(CSV.File(msg)) from the DataFrames and CSV packages, respectively. Alternatively String(msg) will convert it to a string, which you can write to a file.\n\nExample\n\njulia> using PubChemCrawler, CSV, DataFrames\n\njulia> cid = get_cid(name=\"estriol\")\n5756\n\njulia> df = CSV.File(query_substructure(;cid)) |> DataFrame      # on Julia 1.0, use `(;cid=cid)`\n11607×4 DataFrame\n│ Row  │ CID       │ MolecularFormula │ MolecularWeight │ XLogP    │\n│      │ Int64     │ String           │ Float64         │ Float64? │\n├──────┼───────────┼──────────────────┼─────────────────┼──────────┤\n│ 1    │ 5756      │ C18H24O3         │ 288.4           │ 2.5      │\n│ 2    │ 5281904   │ C24H32O9         │ 464.5           │ 1.1      │\n│ 3    │ 27125     │ C18H24O4         │ 304.4           │ 1.5      │\n...\n\nwill query for derivatives of estriol.\n\ninfo: Info\nFor complex queries that risk timing out, consider query_substructure_pug in combination with get_for_cids.\n\n\n\n\n\n","category":"function"},{"location":"#PubChemCrawler.get_for_cids","page":"Home","title":"PubChemCrawler.get_for_cids","text":"msg = get_for_cids(cids; properties|xrefs|cids_type|record_type, output=\"CSV\")\n\nRetrieve the given properties, xrefs, CIDs, or records, respectively, for a list of compounds specified by their cids. The documentation for these traits can be found at http://pubchemdocs.ncbi.nlm.nih.gov/pug-rest; this URL will be referred to as PUGREST below.\n\nproperties include structural features like the molecular formula, number of undefined stereocenters, and so on. Specify these as a comma-separated list from among the choices in PUGREST under \"Compound Property Tables\".\nxrefs (\"cross-references\") include identifiers used by other databases, e.g., the CAS (Registry) number, PubMedID, and so on. The supported values for xrefs are available at PUGREST under \"XRefs\".\ncids_type is used to retrieve CIDs for compounds related to those specified in cids; see PUGREST under \"SIDS / CIDS / AIDS\".\nrecord_type is used to retrieve data files and to specify options for these files, e.g., 2d or 3d SDF files. See PUGREST under \"Full-record Retrieval\".\n\noutput specifies the output format. Not all options are applicable to all queries; for example, \"CSV\" is appropriate for properties queries but \"SDF\" might be used for a record_type query. See PUGREST, \"Output\".\n\nExamples\n\njulia> using PubChemCrawler, CSV, DataFrames, JSON3\n\njulia> cids = [get_cid(name=\"cyclic guanosine monophosphate\"), get_cid(name=\"aspirin\")]\n2-element Array{Int64,1}:\n 135398570\n      2244\n\njulia> CSV.File(get_for_cids(cids; properties=\"MolecularFormula,XLogP\", output=\"CSV\")) |> DataFrame\n2×3 DataFrame\n Row │ CID        MolecularFormula  XLogP\n     │ Int64      String            Float64\n─────┼──────────────────────────────────────\n   1 │ 135398570  C10H12N5O7P          -3.4\n   2 │      2244  C9H8O4                1.2\n\njulia> open(\"/tmp/aspirin_3d.sdf\", \"w\") do io    # save the 3d SDF file for aspirin (CID 2244)\n           write(io, get_for_cids(2244; record_type=\"3d\", output=\"SDF\"))\n       end\n4055\n\njulia> dct = JSON3.read(get_for_cids(cids; xrefs=\"RN,\", output=\"JSON\"));   # get the Registry Number(s) (CAS)\n\njulia> dct[:InformationList][:Information]\n2-element JSON3.Array{JSON3.Object,Array{UInt8,1},SubArray{UInt64,1,Array{UInt64,1},Tuple{UnitRange{Int64}},true}}:\n {\n   \"CID\": 135398570,\n    \"RN\": [\n            \"40732-48-7\",\n            \"7665-99-8\"\n          ]\n}\n {\n   \"CID\": 2244,\n    \"RN\": [\n            \"11126-35-5\",\n            \"156865-15-5\",\n            \"50-78-2\",\n            \"52080-78-1\",\n            \"921943-73-9\",\n            \"98201-60-6\",\n            \"99512-66-0\"\n          ]\n}\n\n\n\n\n\n","category":"function"},{"location":"#Utilities","page":"Home","title":"Utilities","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"parse_formula\natomregex","category":"page"},{"location":"#PubChemCrawler.parse_formula","page":"Home","title":"PubChemCrawler.parse_formula","text":"atomcounts = parse_formula(str::AbstractString)\n\nParse str as a chemical formula, return a list of atom=>multiplicity pairs.\n\nExample\n\njulia> parse_formula(\"C2CaH2O6\")\n4-element Vector{Pair{String, Int64}}:\n  \"C\" => 2\n \"Ca\" => 1\n  \"H\" => 2\n  \"O\" => 6\n\n\n\n\n\n","category":"function"},{"location":"#PubChemCrawler.atomregex","page":"Home","title":"PubChemCrawler.atomregex","text":"rex = atomregex(chemicalsymbol)\n\nCreate a regular expression for detecting how many atoms of type chemicalsymbol are in a molecular formula.\n\nExamples\n\nThe formula for calcium bicarbonate is Ca(HCO3)2, i.e., C2CaH2O6.\n\njulia> match(atomregex(\"C\"), \"C2CaH2O6\")\nRegexMatch(\"C2\", 1=\"2\")\n\njulia> match(atomregex(\"Ca\"), \"C2CaH2O6\")\nRegexMatch(\"Ca\", 1=\"\")\n\njulia> match(atomregex(\"H\"), \"C2CaH2O6\")\nRegexMatch(\"H2\", 1=\"2\")\n\njulia> match(atomregex(\"O\"), \"C2CaH2O6\")\nRegexMatch(\"O6\", 1=\"6\")\n\nNote that the regex for \"C\" does not match \"Ca\".\n\n\n\n\n\n","category":"function"}]
}
