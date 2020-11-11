using PubChemCrawler
using CSV
using JSON3
using DataFrames
using Test

@testset "PubChemCrawler.jl" begin
    # Note that PubChem limits request to no more than 5 per second, so we need to introduce delays in the tests
    cid = get_cid(;name="estriol")
    @test cid == 5756
    sleep(0.2)
    @test get_cid(smiles="CC(=O)OC1=CC=CC=C1C(=O)O") == 2244

    # substructure search
    sleep(0.2)
    df = CSV.File(query_substructure(;cid=cid, output="CSV")) |> DataFrame
    @test !isempty(df)

    # properties
    sleep(0.2)
    df2 = CSV.File(get_for_cids(df.CID[1:2]; properties="MolecularFormula")) |> DataFrame
    @test df2.CID == df.CID[1:2]
    @test parse_formula(df[1,"MolecularFormula"]) == ["C"=>18, "H"=>24, "O"=>3]
    @test parse(Int, match(atomregex("C"), df[1,"MolecularFormula"]).captures[1]) == 18

    # xrefs
    sleep(0.4)  # next one is two requests
    cids = [get_cid(name="cyclic guanosine monophosphate"), get_cid(name="aspirin")]
    sleep(0.2)
    dct = JSON3.read(get_for_cids(cids; xrefs="RN,", output="JSON"))
    @test dct[:InformationList][:Information][1][:RN] == ["40732-48-7", "7665-99-8"]
    sleep(0.2)
    @test chomp(String(get_for_cids(cids[1]; xrefs="RN,", output="TXT"))) == "40732-48-7\n7665-99-8"
end
