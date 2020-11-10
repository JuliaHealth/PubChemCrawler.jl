using PubChemCrawler
using CSV
using DataFrames
using Test

@testset "PubChemCrawler.jl" begin
    cid = get_cid(;name="estriol")
    @test cid == 5756
    sleep(0.1)
    df = CSV.File(query_substructure(;cid, output="CSV")) |> DataFrame
    @test !isempty(df)
    sleep(0.1)
    df2 = CSV.File(get_for_cids(df.CID[1:2])) |> DataFrame
    @test df2.CID == df.CID[1:2]
    @test parse_formula(df[1,"MolecularFormula"]) == ["C"=>18, "H"=>24, "O"=>3]
    @test parse(Int, match(atomregex("C"), df[1,"MolecularFormula"]).captures[1]) == 18
end
