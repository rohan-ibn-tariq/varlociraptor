//! macros.rs
//!
//! Macros for varlociraptor.
//!
//! Sections of Macros included are:
//! 1. String Constants - Macros to define string constants;
//!

/* =============== STRING CONSTANTS ================ */

/// Standard INFO fields always propagated in varlociraptor preprocessing.
///
/// Returns a comma-separated string of field names.
/// These fields are explicitly handled by preprocessing commands.
///
/// # TODO
/// Replace `OMIT_AUX_INFO` in `calling/variants/mod.rs` with this macro.
#[macro_export]
macro_rules! standard_omit_aux_info {
    () => {
        "MATEID, EVENT, SVLEN, SVTYPE, END"
    };
}

/// MSI-specific INFO fields always propagated in MSI preprocessing,
/// extending the standard fields with MSI-required fields.
#[macro_export]
macro_rules! preprocess_msi_omit_aux_info {
    () => {
        concat!(standard_omit_aux_info!(), ", HETEROZYGOSITY")
    };
}

/* ================================================= */
