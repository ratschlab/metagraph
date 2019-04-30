//
//  tests.hpp
//  lossless_dbg
//
//  Created by Jan Studený on 11/03/2019.
//  Copyright © 2019 Jan Studený. All rights reserved.
//

#ifndef tests_h
#define tests_h

#include <gtest/gtest.h>
#define private public
#define protected public
#include "utilities.hpp"
#include "path_database_list_of_bifurcation_choices.hpp"
#include "samplers.hpp"
#include "path_database_baseline.hpp"
#include "path_database_baseline_wavelet_deprecated.hpp"
#include "path_database_baseline_wavelet.hpp"
#include <filesystem>
using namespace std;
namespace fs = std::filesystem;
//#include "path_database.hpp"


const int test_seed = 3424;
const string random_generator_state = "3424 3838682849 1000327468 417786463 3810545279 1529041649 1625434294 2449801018 3809847584 2466543064 2496276748 3944525969 2295188390 3702004673 2675688856 2278186961 2196027727 2558677106 1535487042 2663502210 370485140 2375209337 191133853 3413809928 3315928687 1137343669 1955301150 1815406166 8679279 4092288488 3867922645 3318082189 1132300326 4218721924 2873008741 4184879558 3964211165 1768100027 2233410440 1611915161 3349647136 514470904 985544962 1620873205 3329908592 3951113100 3413761689 2372615409 1569855247 2723911351 2361152155 481261456 1025915908 3931099593 491314920 2237798847 2409974985 3755798352 495113977 2524516728 443821918 2428632659 232369971 1141467742 3728527035 2141587161 2249345914 1608711835 1023728902 1079863715 2647678768 2022002945 468841800 2544979377 1601008361 1023704531 4007694219 3089131765 263624897 2215677300 2624976350 1679866909 4200149854 624494852 611300072 3945250269 3831215852 27136930 491614530 4059562339 3107253050 1917437299 266528854 1429581899 1670129552 2288989588 380070542 1509230695 3790923424 1199651250 1975619331 3421732399 3368702914 887903372 3125823396 1972247271 834646824 4254734131 1719916124 407263262 1236980292 3978224296 1797910759 1472274735 1245033880 2291042000 1847563086 3010684320 2908277600 3498946593 1142337506 2487803912 3639328364 3337914694 953357749 1440247270 646249377 713519108 41768468 3641682533 3056128192 1126182669 238287936 1379976645 901769946 1297896073 2404341808 2178846275 2796939567 2076775500 2085435373 1960336553 1267365590 3675976290 1778289109 1039848757 3766647675 1742073579 3030940902 193677961 4002855331 1580063927 2217704038 402050573 2801746107 1298565784 3218051065 2496536996 4175371548 1987384538 1515679239 4169768703 1564775950 4005968782 1726936133 3671769977 994142664 2799879311 243401033 553227126 1654346808 1625930792 3272093145 2344206767 3500801775 3572621259 4071051928 3983903448 3333159961 1215979765 2899760376 1938936151 2733153700 430836021 3775528865 2016251299 36869028 678626415 1976117127 3140308635 3391022107 2733168183 401727657 1413944430 2902685581 3057824558 275201312 3830529637 1262236932 1984141504 1857129965 578046693 2084117027 2691793461 3472412543 3415685305 1159392304 3445385508 3452722227 1413877697 743569298 3067203437 3906092703 3837276257 2034076032 4066658748 3846262835 3623769033 3568576908 3426413894 3038546453 885639920 2647097742 2665351451 1277622461 60866317 704927235 2211613458 1005350452 514674537 1371440723 2051727681 2060853032 860241174 3305169048 3116929298 1496501308 2986490366 1372372314 1754884310 4018302147 1139255473 2692464226 1045658835 2994875187 4040115274 2538603203 55664668 1280174084 1074411506 3349028569 802871037 2605660109 2244463784 83933712 860617039 177054763 2548680952 44045732 2084130487 2629137618 3771361045 446657204 3895891467 1309375024 4064440158 1767486907 45825389 4055470605 2478972563 803776835 1312072574 314341675 3043021320 1440758020 2333945100 2734769562 82234893 2444181815 1753028608 1060766333 3399672170 2462261383 1392902804 3872637669 3240295387 4001101142 4270147752 3743821719 915602053 632358043 512410954 3982127446 11023534 525330380 2217124003 530221485 3948996458 1492221591 304035929 4058504521 3404523871 165564026 2047492433 741316032 4264154609 2610817964 3257932505 256852790 330787715 612075237 3425683088 2451768375 3733194530 2412673599 1334566988 3679696029 3959476627 124305166 468967365 3294279417 2089548771 1211388780 3242541124 94985287 1478781000 3895146259 4210313879 945201580 2456213029 440294829 1182239628 330532077 1987135438 2764578041 3377644374 3951202009 2833975123 3403365703 3106181415 823988717 1753934038 2806639913 809735758 977013790 419354415 4093382117 2401997593 933723651 1301894028 1880059135 2589716373 1520326387 4120760539 837746330 2277563941 1292600775 3230449283 2130291174 3451758730 314822773 1090706578 976269929 831493592 2646521764 1511017451 208940736 1578840879 2103496854 2935233540 4118300112 580286642 2429739438 3115062865 4179376949 3394061253 2351608982 2273866717 3186691957 2443364974 2894178840 300287935 4095690457 4114063233 3084681674 2011905385 3246646154 4069267344 1390934403 2224686031 1399273063 2828996037 4187025419 4141546673 2524902852 2449993385 2257038851 3719114994 656165283 1493584350 2518449547 3005929886 1637898270 3494141134 3325907829 404100387 2851804005 1765918778 414749663 2802691988 1350012360 1720353512 1964202121 1455077189 3600397170 2143438388 2796229001 1823363192 1764065887 3340943033 3621571334 4226074526 2170603415 655391856 3811117528 1047795728 756944890 3409910093 2375138162 937288925 1209045471 4054539333 451338062 767089015 479165605 315849676 222561584 116178341 946282703 2867313250 147450264 148079281 3946465935 431511287 1627750447 570831331 3730385741 2289882501 1916174083 2324049803 2401244111 2832959140 3493132098 2885572458 1006052814 3697224205 3983480654 1063741482 3954218076 2661966406 665528992 3937891053 1399044020 2778356024 3627761586 1850530726 460327349 1750643516 1167416805 610865353 3639489059 3987531895 2563361692 2129061167 1308625664 2230248000 1361221606 557944576 493120734 4176724853 3422569838 2201860834 1820105794 111051602 2685317950 979096721 653490203 1444009870 3125711955 3365330654 886387739 3250389906 3681682209 93745239 2094654273 324266543 3283586427 4033350217 110340900 2094320423 3065347570 3884982693 1902714484 821601056 1464106392 275476566 1406884328 973548520 3467345028 2868732480 2539329544 1460918257 552838064 478725745 2369766807 502837964 3730946176 1474357428 3241651823 1547935395 1157646834 3931453928 3372981953 3101259413 1080269471 2539425891 527425875 1946047182 1968207035 4954483 1226250097 364327235 930630787 2667400388 3916382772 2799014090 1916086528 217167998 2861251024 4154274805 93202474 364532911 1855853097 2603389927 3137267833 1884851624 3797833423 994012319 1568666591 2953774779 3147147811 2902379052 2314888270 1032470565 2699173315 3077113936 1679961222 4144026224 689411213 1621627600 3289925541 1500912815 2981392600 3129684533 2857058535 3618515342 3939359703 1084764379 3530945082 3535279030 2334041507 3173090240 643014 2054811995 633038016 3977871871 260318892 1988710429 1087347022 807126382 2778645674 2407104653 2024096177 355081399 3757852027 3858539169 3122847284 4096168857 4184582158 415450222 2343231412 1914481437 3668279388 1805090508 2033958451 414764557 3314106997 4243742947 2811763382 809033563 4118074431 2052053509 2739481070 1220597879 2602995434 2701691877 1653527425 436794335 3475214427 3551464729 2244548516 327691489 3135578153 1607054940 4122458647 3399295051 1211909840 3714370526 417904027 3090979730 538538556 3083874329 2745069589 471829122";

TEST(SamplerTest,SampleNoRandom) {
    auto state = stringstream(random_generator_state);
    auto generator = mt19937();
    state >> generator;
    auto sampler = NoisySampler("AAAAAAAAA",generator);
    ASSERT_EQ(sampler.sample(2),"AA");
}

TEST(SamplerTest,DISABLED_SampleNormal) {
    // TODO: figure out how to make random generator platform independent
    auto state = stringstream(random_generator_state);
    auto generator = mt19937();
    state >> generator;
    auto sampler = NoisySampler("ADFAGADFDS",generator);
    ASSERT_EQ(sampler.sample(4),"ADFD");
}
TEST(SamplerTest,SampleCoverage) {
    auto state = stringstream(random_generator_state);
    auto generator = mt19937();
    state >> generator;
    auto sequence = "ADFAGADFDS"s;
    auto sampler = NoisySampler(sequence,generator);
    auto reads = sampler.sample_coverage(sequence.length()/2, 1);
    ASSERT_EQ(reads.size(), 2);
}
TEST(SamplerTest,SubSample) {
    auto state = stringstream(random_generator_state);
    auto generator = mt19937();
    state >> generator;
    auto sequence = "ADFAGADFDSAFDAGDGDASGDSFADSGDSVBABDFDASFDSAFADAFDSFASDVCZVCXFDAF"s;
    auto sampler = SubSampler(sequence,10,generator);
    auto reads = sampler.sample_coverage(5, 1);
    ASSERT_EQ(reads.size(), 2);
}

// Depends on large file -> tested and works
//TEST(CompressingReads,GetChromosomeWorks) {
//    auto chromosome = get_human_chromosome(CHROMOSOME_NUMBER);
//    EXPECT_EQ(chromosome.length(), 133'797'422);
//    EXPECT_EQ(chromosome.substr(0,10),"NNNNNNNNNN");
//}


vector<string> reads_for_testing_short = {"ATGCGATCGATATGCGAGA",
                                          "ATGCGATCGAGACTACGAG",
                                          "GTACGATAGACATGACGAG",
                                          "ACTGACGAGACACAGATGC"};

template <typename T>
void check_compression_decompression(T& db,vector<string>& reads) {
    auto handles = db.encode(reads);
    vector<string> decompressed_reads;
    decompressed_reads.reserve(handles.size());
    for(auto& handle : handles) {
        decompressed_reads.push_back(db.decode(handle));
    }
    ASSERT_EQ(reads,decompressed_reads);
}

template <typename T>
void check_decompression_all(T& db,vector<string>& reads) {
    db.encode(reads);
    auto decompressed_reads = db.decode_all_reads();
    ASSERT_EQ(multiset<string>(all(reads)),multiset<string>(all(decompressed_reads)));
}

template <typename T>
void check_decompression_inverse(T& db,vector<string>& reads) {
    db.encode(reads);
    auto decompressed_reads = db.decode_all_reads_inverse();
    ASSERT_EQ(multiset<string>(all(reads)),multiset<string>(all(decompressed_reads)));
}

template <typename T>
void check_compression_decompression(vector<string>& reads, int k_kmer=21) {
    auto db = T(reads,k_kmer);
    check_compression_decompression(db,reads);
}

template <typename T>
void short_reads_decode_all() {
    auto db = T(reads_for_testing_short,5);
    check_decompression_all(db,reads_for_testing_short);
}

template <typename T>
void short_reads_decode_inverse() {
    auto db = T(reads_for_testing_short,5);
    check_decompression_inverse(db,reads_for_testing_short);
}

template <typename T>
void short_identity_test() {
    check_compression_decompression<T>(reads_for_testing_short,5);
}

template <class T=PathDatabaseBaselineWavelet<>>
void check_paths_going_through() {
    auto db = T(reads_for_testing_short,5);
    db.encode(reads_for_testing_short);
    ASSERT_EQ(db.get_paths_going_through(db.graph.kmer_to_node("CATGA")),vector<typename T::path_id>({{db.graph.kmer_to_node("GTACG"),0}}));
}
template<typename T=PathDatabaseBaselineWavelet<>>
void check_small_get_next_consistent_node() {
    string middle = "ACTGCGT";
    vector<string> reads = {"C" + middle + "T","A" + middle + "G"};
    auto db = T(reads,5);
    db.encode(reads);
    auto split_node = db.graph.kmer_to_node("TGCGT");
    ASSERT_EQ(db.get_next_consistent_node(split_node,"C" + middle),db.graph.kmer_to_node("GCGT"s + "T"s));
    ASSERT_EQ(db.get_next_consistent_node(split_node,"A" + middle),db.graph.kmer_to_node("GCGT"s + "G"s));
    ASSERT_EQ(db.get_next_consistent_node(split_node,middle),0);
}


template <typename T>
void serialization_deserialization_test(vector<string>& reads, int k_kmer=21) {
    auto db = T(reads,k_kmer);
    auto handles = db.encode(reads);
    vector<string> decompressed_reads;
    decompressed_reads.reserve(handles.size());
    auto output_folder = fs::temp_directory_path() / "serdes/" ;
    cout << "Files saved to: " << output_folder << endl;
    fs::create_directories(output_folder);
    db.serialize(output_folder);
    auto newdb = T::deserialize(output_folder);
    decompressed_reads.reserve(handles.size());
    for(auto& handle : handles) {
        decompressed_reads.push_back(newdb.decode(handle));
    }
    ASSERT_EQ(reads,decompressed_reads);
}

template <typename T>
void short_serdes_test() {
    serialization_deserialization_test<T>(reads_for_testing_short,5);
}

TEST(PathDatabase,IncomingTable) {
    DBGSuccinct graph = DBGSuccinct(21);
    IncomingTable table(graph);
    table.joins = IncomingTable<>::bit_vector_t({1,0,0,1,0,1,1,0,1});
    vector<int> v = {1,2,3,4};
    table.edge_multiplicity_table = sdsl::enc_vector<>(v);
    ASSERT_EQ(table.branch_size_rank(1,0),1);
    ASSERT_EQ(table.branch_size_rank(1,1),2);
    ASSERT_EQ(table.branch_size_rank(2,0),3);
    ASSERT_EQ(table.size(3),0);
    ASSERT_EQ(table.branch_size_rank(4,0),4);
}

TEST(PathDatabase,IdentityTestCompressedReads) {
    short_identity_test<PathDatabaseListBC>();
}

TEST(PathDatabase,IdentityTestPathDatabaseBaseline) {
    short_identity_test<PathDatabaseBaseline<>>();
}

TEST(PathDatabase,IdentityTestPathDatabaseBaselineWaveletDeprecated) {
    short_identity_test<PathDatabaseBaselineWaveletDeprecated>();
}

TEST(PathDatabase,IdentityTestPathDatabaseBaselineWavelet) {
    short_identity_test<PathDatabaseBaselineWavelet<>>();
}

TEST(PathDatabase,SerDesTest) {
    short_serdes_test<PathDatabaseBaselineWaveletDeprecated>();
}

TEST(PathDatabase,DecodeAllInverse) {
    short_reads_decode_inverse<PathDatabaseBaselineWavelet<>>();
}

TEST(PathDatabase,DecodeAll) {
    short_reads_decode_all<PathDatabaseBaselineWavelet<>>();
}
TEST(PathDatabase,PathsGoingThrough) {
    check_paths_going_through<PathDatabaseBaselineWavelet<>>();
}
TEST(PathDatabase,ConsistentNode) {
    check_small_get_next_consistent_node<>();
}

#if defined(__linux__) || false

template <typename T>
void long_identity_test() {
    string reads_filename = "/cluster/home/studenyj/genomic_data/human_chr10_artifical_reads.fasta";
    auto reads = read_reads_from_fasta(reads_filename);
    check_compression_decompression<T>(reads);
}
TEST(PathDatabase,LongTestCompressedReads) {
    long_identity_test<PathDatabaseListBC>();
}

TEST(PathDatabase,LongTestBaseline) {
    long_identity_test<PathDatabaseBaseline<>>();
}

#endif

#endif /* tests_h */
