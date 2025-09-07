using PubChemCrawler
using CSV
using JSON3
using DataFrames
using BrokenRecord: BrokenRecord, playback
using HTTP    # needed to make BSON happy upon playback
using Test

## NOTE: the interactions with PubChem are recorded using BrokenRecord and stored in the test/http_record
## directory. If you want to check that this still works with the "real" PubChem server, just delete
## that directory and re-run the tests.

const allrecordings = [joinpath("http_record", file) for file in [
    "estriol_cid.bson",
    "aspirin_cid_from_smiles.bson",
    "estriol_substructure.bson",
    "smarts.bson",
    "smarts_pug.bson",
    "aspirin_smiles_from_cid.bson",
    "estriol_substructure_formulas.bson",
    "asprin_sdf_3d.bson",
    "cGMP_cid.bson",
    "aspirin_cid_from_name.bson",
    "CAS_as_json.bson",
    "CAS_as_txt.bson",
    "estriol_synonyms_from_name.bson",
    "estriol_synonyms_from_cid.bson",
    "estriol_synonyms_from_smiles.bson"
                           
]]
const get_recordings = !all(isfile, allrecordings)

if !isdir("http_record")
    mkdir("http_record")
end
BrokenRecord.configure!(; path="http_record")

@testset "PubChemCrawler.jl" begin
    # Note that PubChem limits requests to no more than 5 per second, so we need to introduce delays in the tests.
    # Moreover, multiple tests being run simultaneously by GitHubActions may count as originating from the same place,
    # so let's be quite conservative about the delays.
    cid_estriol = playback(() -> get_cid(;name="estriol"), "estriol_cid.bson")
    @test cid_estriol == 5756
    sleep(2.0 * get_recordings)
    @test playback(() -> get_cid(smiles="CC(=O)OC1=CC=CC=C1C(=O)O"), "aspirin_cid_from_smiles.bson") == 2244
    cids_aspirin = playback(() -> get_cids(;cas_number="50-78-2"), "aspirin_cids.bson")
    @test cids_aspirin == [2244, 67252, 3434975, 12280114]

    # substructure search
    sleep(2.0 * get_recordings)
    df = CSV.File(playback(() -> query_substructure(;cid=cid_estriol, output="CSV"), "estriol_substructure.bson")) |> DataFrame
    @test 27125 ∈ df.CID   # check that estetrol has estriol as a substructure
    sleep(2.0 * get_recordings)
    df13 = CSV.File(playback(() -> query_substructure(;smarts="[r13]Br", output="CSV"), "smarts.bson")) |> DataFrame  # brominated 13-atom ring structures
    @test 118303825 ∈ df13.CID
    # The recommended approach for substructure searches is `query_substructure_pug`
    sleep(5.0 * get_recordings)
    cids13 = playback(() -> query_substructure_pug(;smarts="[r13]Br", poll_interval=10*get_recordings), "smarts_pug.bson")  # brominated 13-atom ring structures, via PUG interface
    @test 118303825 ∈ cids13

    # properties
    sleep(5.0 * get_recordings)
    smiles = String(playback(() -> get_for_cids(2244, properties="CanonicalSMILES", output="TXT"), "aspirin_smiles_from_cid.bson"))
    @test chomp(smiles) == "CC(=O)OC1=CC=CC=C1C(=O)O"
    sleep(5.0 * get_recordings)
    df2 = CSV.File(playback(() -> get_for_cids(df.CID[1:2]; properties="MolecularFormula"), "estriol_substructure_formulas.bson")) |> DataFrame
    @test df2.CID == df.CID[1:2]
    @test parse_formula(df[1,"MolecularFormula"]) == ["C"=>18, "H"=>24, "O"=>3]
    @test parse(Int, match(atomregex("C"), df[1,"MolecularFormula"]).captures[1]) == 18

    # structure files (SDF)
    sleep(5.0 * get_recordings)
    sdf = String(playback(() -> get_for_cids(2244, output="SDF", record_type="3d"), "asprin_sdf_3d.bson"))
    @test occursin("PUBCHEM_MMFF94_PARTIAL_CHARGES", sdf)
    line = split(sdf, '\n')[5]
    flds = split(line)
    @test parse(Float32, flds[1]) != 0 && parse(Float32, flds[2]) != 0 && parse(Float32, flds[3]) != 0

    # parent compounds (sodium acetate is 517045, acetic acid is 176)
    sleep(5.0 * get_recordings)
    str = String(playback(() -> get_for_cids(517045; cids_type="parent", output="TXT"), "sodium_acetate_parent.bson"))
    @test parse(Int, chomp(str)) == 176

    # xrefs
    sleep(4.0 * get_recordings)  # next one is two requests
    cids = [playback(() -> get_cid(name="cyclic guanosine monophosphate"), "cGMP_cid.bson")
            playback(() -> get_cid(name="aspirin"), "aspirin_cid_from_name.bson")]
    sleep(5.0 * get_recordings)
    dct = JSON3.read(playback(() -> get_for_cids(cids; xrefs="RN,", output="JSON"), "CAS_as_json.bson"))
    @test dct[:InformationList][:Information][1][:RN] ==  ["231-641-6", "40732-48-7", "7665-99-8"]
    sleep(5.0 * get_recordings)
    @test chomp(String(playback(() -> get_for_cids(cids[1]; xrefs="RN,", output="TXT"), "CAS_as_txt.bson"))) == "231-641-6\n40732-48-7\n7665-99-8"

    # pug
    cid =  playback(() -> pug(:compound, :name, "ethanol", :cids, :txt, return_text = true), "ethanol_pug.bson")
    @test cid == "702"

    #synonyms
    @test "Trimesta" ∈ playback(() -> get_synonyms(name="estriol"), "estriol_synonyms_from_name.bson")
    sleep(5.0 * get_recordings)
    @test "Trimesta" ∈ playback(() -> get_synonyms(cid=cid_estriol), "estriol_synonyms_from_cid.bson")
    sleep(5.0 * get_recordings)
    @test "Trimesta" ∈ playback(() -> get_synonyms(smiles="C[C@]12CC[C@H]3[C@H]([C@@H]1C[C@H]([C@@H]2O)O)CCC4=C3C=CC(=C4)O"), "estriol_synonyms_from_smiles.bson")
end
