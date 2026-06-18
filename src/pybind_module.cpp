#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "aligner.h"
#include "sequence.h"
#include "fasta_reader.h"
#include "fastq_reader.h"
#include "fm_index.h"
#include "suffix_array.h"
#include "bwt.h"
#include "backward_search.h"
#include "sam_writer.h"

namespace py = pybind11;
using namespace bwt_aligner;

PYBIND11_MODULE(bwt_aligner, m) {
    m.doc() = "High-performance BWT-based genomic sequence aligner";

    py::class_<ReferenceSequence>(m, "ReferenceSequence")
        .def(py::init<>())
        .def_readwrite("name", &ReferenceSequence::name)
        .def_readwrite("sequence", &ReferenceSequence::sequence)
        .def_readwrite("length", &ReferenceSequence::length)
        .def_readwrite("global_offset", &ReferenceSequence::global_offset);

    py::class_<FastqRead>(m, "FastqRead")
        .def(py::init<>())
        .def_readwrite("id", &FastqRead::id)
        .def_readwrite("sequence", &FastqRead::sequence)
        .def_readwrite("quality", &FastqRead::quality)
        .def_readwrite("is_first", &FastqRead::is_first)
        .def_readwrite("pair_id", &FastqRead::pair_id);

    py::class_<AlignmentResult>(m, "AlignmentResult")
        .def(py::init<>())
        .def_readwrite("read_id", &AlignmentResult::read_id)
        .def_readwrite("ref_name", &AlignmentResult::ref_name)
        .def_readwrite("ref_id", &AlignmentResult::ref_id)
        .def_readwrite("position", &AlignmentResult::position)
        .def_readwrite("is_reverse", &AlignmentResult::is_reverse)
        .def_readwrite("mapq", &AlignmentResult::mapq)
        .def_readwrite("cigar", &AlignmentResult::cigar)
        .def_readwrite("sequence", &AlignmentResult::sequence)
        .def_readwrite("quality", &AlignmentResult::quality)
        .def_readwrite("is_secondary", &AlignmentResult::is_secondary)
        .def_readwrite("is_supplementary", &AlignmentResult::is_supplementary)
        .def_readwrite("is_paired", &AlignmentResult::is_paired)
        .def_readwrite("is_proper_pair", &AlignmentResult::is_proper_pair)
        .def_readwrite("mate_is_reverse", &AlignmentResult::mate_is_reverse)
        .def_readwrite("mate_ref_id", &AlignmentResult::mate_ref_id)
        .def_readwrite("mate_position", &AlignmentResult::mate_position)
        .def_readwrite("insert_size", &AlignmentResult::insert_size)
        .def_readwrite("is_unmapped", &AlignmentResult::is_unmapped)
        .def_readwrite("mate_is_unmapped", &AlignmentResult::mate_is_unmapped)
        .def_readwrite("num_mismatches", &AlignmentResult::num_mismatches);

    py::class_<SearchRange>(m, "SearchRange")
        .def(py::init<>())
        .def_readwrite("low", &SearchRange::low)
        .def_readwrite("high", &SearchRange::high)
        .def("valid", &SearchRange::valid)
        .def("count", &SearchRange::count);

    py::class_<SeedHit>(m, "SeedHit")
        .def(py::init<>())
        .def_readwrite("ref_position", &SeedHit::ref_position)
        .def_readwrite("query_start", &SeedHit::query_start)
        .def_readwrite("query_end", &SeedHit::query_end)
        .def_readwrite("num_mismatches", &SeedHit::num_mismatches)
        .def_readwrite("is_reverse", &SeedHit::is_reverse);

    py::class_<AlignerConfig>(m, "AlignerConfig")
        .def(py::init<>())
        .def_readwrite("num_threads", &AlignerConfig::num_threads)
        .def_readwrite("seed_length", &AlignerConfig::seed_length)
        .def_readwrite("max_mismatches", &AlignerConfig::max_mismatches)
        .def_readwrite("min_seed_hits", &AlignerConfig::min_seed_hits)
        .def_readwrite("batch_size", &AlignerConfig::batch_size)
        .def_readwrite("min_mapq", &AlignerConfig::min_mapq);

    py::class_<FastaReader>(m, "FastaReader")
        .def(py::init<>())
        .def("open", &FastaReader::open)
        .def("close", &FastaReader::close)
        .def("read_all", &FastaReader::read_all)
        .def("total_length", &FastaReader::total_length)
        .def("num_sequences", &FastaReader::num_sequences);

    py::class_<FastqReader>(m, "FastqReader")
        .def(py::init<>())
        .def("open", &FastqReader::open)
        .def("close", &FastqReader::close)
        .def("read_next", &FastqReader::read_next)
        .def("read_batch", &FastqReader::read_batch)
        .def("total_reads", &FastqReader::total_reads)
        .def("is_open", &FastqReader::is_open);

    py::class_<PairedFastqReader>(m, "PairedFastqReader")
        .def(py::init<>())
        .def("open", &PairedFastqReader::open)
        .def("close", &PairedFastqReader::close)
        .def("read_batch_paired", &PairedFastqReader::read_batch_paired)
        .def("total_reads", &PairedFastqReader::total_reads);

    py::class_<SuffixArray>(m, "SuffixArray")
        .def(py::init<>())
        .def("build", &SuffixArray::build)
        .def("get_sa", &SuffixArray::get_sa, py::return_value_policy::reference_internal)
        .def("size", &SuffixArray::size)
        .def("clear", &SuffixArray::clear)
        .def("__getitem__", [](const SuffixArray& self, size_t i) {
            if (i >= self.size()) throw py::index_error();
            return self[i];
        });

    py::class_<BWT>(m, "BWT")
        .def(py::init<>())
        .def("build", &BWT::build)
        .def("get_bwt", &BWT::get_bwt, py::return_value_policy::reference_internal)
        .def("size", &BWT::size)
        .def("get_primary", &BWT::get_primary)
        .def("get_primary_pos", &BWT::get_primary_pos)
        .def("clear", &BWT::clear)
        .def("__getitem__", [](const BWT& self, size_t i) {
            if (i >= self.size()) throw py::index_error();
            return self[i];
        });

    py::class_<FMIndex>(m, "FMIndex")
        .def(py::init<>())
        .def("build", py::overload_cast<const std::string&>(&FMIndex::build))
        .def("get_c", &FMIndex::get_c)
        .def("get_occ", &FMIndex::get_occ)
        .def("get_sa", &FMIndex::get_sa, py::return_value_policy::reference_internal)
        .def("get_bwt", &FMIndex::get_bwt, py::return_value_policy::reference_internal)
        .def("get_text", &FMIndex::get_text, py::return_value_policy::reference_internal)
        .def("size", &FMIndex::size)
        .def("get_sa_value", &FMIndex::get_sa_value)
        .def("resolve_position", &FMIndex::resolve_position)
        .def("clear", &FMIndex::clear);

    py::class_<BackwardSearch>(m, "BackwardSearch")
        .def(py::init<>())
        .def("set_fm_index", &BackwardSearch::set_fm_index)
        .def("search", &BackwardSearch::search)
        .def("search_reverse", &BackwardSearch::search_reverse)
        .def("find_seeds", &BackwardSearch::find_seeds)
        .def("find_seeds_both_strands", &BackwardSearch::find_seeds_both_strands)
        .def("verify_alignment", [](const BackwardSearch& self,
                                     const std::string& query,
                                     int64_t ref_position,
                                     bool is_reverse,
                                     int32_t max_mismatches) {
            int32_t mm = 0;
            bool ok = self.verify_alignment(query, ref_position, is_reverse, max_mismatches, mm);
            return py::make_tuple(ok, mm);
        });

    py::class_<SamWriter>(m, "SamWriter")
        .def(py::init<>())
        .def("open", &SamWriter::open)
        .def("close", &SamWriter::close)
        .def("write_header", &SamWriter::write_header,
             py::arg("references"),
             py::arg("program_id") = "bwt_aligner",
             py::arg("program_version") = "1.0.0")
        .def("write_alignment", &SamWriter::write_alignment)
        .def("write_alignments", &SamWriter::write_alignments)
        .def("total_written", &SamWriter::total_written);

    py::class_<Aligner>(m, "Aligner")
        .def(py::init<>())
        .def("set_config", &Aligner::set_config)
        .def("get_config", &Aligner::get_config, py::return_value_policy::reference_internal)
        .def("build_index", &Aligner::build_index)
        .def("load_reference", &Aligner::load_reference)
        .def("align_single", &Aligner::align_single)
        .def("align_paired", &Aligner::align_paired)
        .def("get_references", &Aligner::get_references, py::return_value_policy::reference_internal)
        .def("get_fm_index", &Aligner::get_fm_index, py::return_value_policy::reference_internal)
        .def("align_read", &Aligner::align_read);

    m.def("reverse_complement", &reverse_complement);
    m.def("complement_base", &complement_base);
    m.def("base_to_index", &base_to_index);
    m.def("index_to_base", &index_to_base);
}
