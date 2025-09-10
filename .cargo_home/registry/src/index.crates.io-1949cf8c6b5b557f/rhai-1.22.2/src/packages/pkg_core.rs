#[cfg(feature = "no_std")]
use std::prelude::v1::*;

use super::*;
use crate::def_package;

def_package! {
    /// Core package containing basic facilities.
    ///
    /// # Contents
    ///
    /// * [`LanguageCorePackage`][super::LanguageCorePackage]
    /// * [`ArithmeticPackage`][super::ArithmeticPackage]
    /// * [`BasicStringPackage`][super::BasicStringPackage]
    /// * [`BasicIteratorPackage`][super::BasicIteratorPackage]
    /// * [`BasicFnPackage`][super::BasicFnPackage]
    /// * [`DebuggingPackage`][super::DebuggingPackage]
    pub CorePackage(lib) :
            LanguageCorePackage,
            ArithmeticPackage,
            BasicStringPackage,
            BasicIteratorPackage,
            BasicFnPackage,
            #[cfg(feature = "debugging")] DebuggingPackage
        {
        lib.set_standard_lib(true);
    }
}
