using PubChemCrawler
using CSV
using JSON3
using DataFrames
using Test

const isci = parse(Bool, get(ENV, "CI", "false"))

@testset "PubChemCrawler.jl" begin
    # Note that PubChem limits requests to no more than 5 per second, so we need to introduce delays in the tests.
    # Moreover, multiple tests being run simultaneously by GitHubActions may count as originating from the same place,
    # so let's be quite conservative about the delays.
    cid = get_cid(;name="estriol")
    @test cid == 5756
    sleep(2.0)
    @test get_cid(smiles="CC(=O)OC1=CC=CC=C1C(=O)O"; verbose=isci*3) == 2244

    # substructure search
    sleep(2.0)
    df = CSV.File(query_substructure(;cid=cid, output="CSV", verbose=isci*3)) |> DataFrame
    @test !isempty(df)

    # properties
    sleep(2.0)
    df2 = CSV.File(get_for_cids(df.CID[1:2]; properties="MolecularFormula", verbose=isci*3)) |> DataFrame
    @test df2.CID == df.CID[1:2]
    @test parse_formula(df[1,"MolecularFormula"]) == ["C"=>18, "H"=>24, "O"=>3]
    @test parse(Int, match(atomregex("C"), df[1,"MolecularFormula"]).captures[1]) == 18

    # xrefs
    sleep(4.0)  # next one is two requests
    cids = [get_cid(name="cyclic guanosine monophosphate", verbose=isci*3), get_cid(name="aspirin", verbose=isci*3)]
    sleep(5.0)
    dct = JSON3.read(get_for_cids(cids; xrefs="RN,", output="JSON", verbose=isci*3))
    @test dct[:InformationList][:Information][1][:RN] == ["40732-48-7", "7665-99-8"]
    sleep(5.0)
    @test chomp(String(get_for_cids(cids[1]; xrefs="RN,", output="TXT", verbose=isci*3))) == "40732-48-7\n7665-99-8"
end
