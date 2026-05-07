//! intersection.rs
//!
//! Streaming intersection alongside with dummy indel
//! injection and variant annotation for MSI detection.
//!
//! This module provides:
//! 1. Streaming intersection of BED regions with VCF variants;
//! 2. Perfect microsatellite repeat detection;
//! 3. Variant filtering;
//! 4. Dummy indel injection for MS regions without variants;
//! 5. MS region annotation (INFO/MS_REGION).
//!
//! The streaming approach processes BED regions sequentially while maintaining
//! a sliding window of VCF variants, enabling memory-efficient analysis of
//! large datasets.
//!
//! Note: This module assumes that both the BED and VCF files are sorted by chromosome
//! (lexicographically) and position, which is a common requirement for genomic analyses.

use std::collections::VecDeque;
use std::path::Path;

use anyhow::{Context, Result};
use bio::io::bed;
use log::{debug, info};
use rust_htslib::bcf::{self, header::HeaderView, Read, Writer};

use crate::errors::Error;
use crate::utils::bcf_utils::{
    get_chrom, get_svlen, is_breakend, is_reference_allele, is_spanning_deletion, is_symbolic,
};
use crate::utils::genomics::{
    calculate_anchor_length, calculate_indel_position, is_clean_indel, is_indel,
};
use crate::utils::ms_bed::{parse_bed_record, BedRegion};
