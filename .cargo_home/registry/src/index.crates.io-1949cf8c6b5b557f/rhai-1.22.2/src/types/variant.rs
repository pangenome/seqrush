//! [`Variant`] trait to to allow custom type handling.

use crate::func::SendSync;
use std::any::{type_name, Any, TypeId};
#[cfg(feature = "no_std")]
use std::prelude::v1::*;

mod private {
    use crate::func::SendSync;
    use std::any::Any;

    /// A sealed trait that prevents other crates from implementing [`Variant`][super::Variant].
    pub trait Sealed {}

    impl<T: Any + Clone + SendSync> Sealed for T {}
}

/// _(internals)_ Trait to represent any type.
/// Exported under the `internals` feature only.
///
/// This trait is sealed and cannot be implemented.
///
/// Currently, [`Variant`] is not [`Send`] nor [`Sync`], so it can practically be any type.
/// Turn on the `sync` feature to restrict it to only types that implement [`Send`] `+` [`Sync`].
#[cfg(not(feature = "sync"))]
pub trait Variant: Any + private::Sealed {
    /// Convert this [`Variant`] trait object to [`&dyn Any`][Any].
    #[must_use]
    fn as_any(&self) -> &dyn Any;

    /// Convert this [`Variant`] trait object to [`&mut dyn Any`][Any].
    #[must_use]
    fn as_any_mut(&mut self) -> &mut dyn Any;

    /// Convert this [`Variant`] trait object to [`Box<dyn Any>`][Any].
    #[must_use]
    fn as_boxed_any(self: Box<Self>) -> Box<dyn Any>;

    /// Get the name of this type.
    #[must_use]
    fn type_name(&self) -> &'static str;

    /// Clone this [`Variant`] trait object.
    #[must_use]
    fn clone_object(&self) -> Box<dyn Variant>;
}

/// _(internals)_ Trait to represent any type.
/// Exported under the `internals` feature only.
///
/// This trait is sealed and cannot be implemented.
#[cfg(feature = "sync")]
pub trait Variant: Any + Send + Sync + private::Sealed {
    /// Convert this [`Variant`] trait object to [`&dyn Any`][Any].
    #[must_use]
    fn as_any(&self) -> &dyn Any;

    /// Convert this [`Variant`] trait object to [`&mut dyn Any`][Any].
    #[must_use]
    fn as_any_mut(&mut self) -> &mut dyn Any;

    /// Convert this [`Variant`] trait object to [`Box<dyn Any>`][Any].
    #[must_use]
    fn as_boxed_any(self: Box<Self>) -> Box<dyn Any>;

    /// Get the name of this type.
    #[must_use]
    fn type_name(&self) -> &'static str;

    /// Clone this [`Variant`] trait object.
    #[must_use]
    fn clone_object(&self) -> Box<dyn Variant>;
}

impl<T: Any + Clone + SendSync> Variant for T {
    #[inline(always)]
    fn as_any(&self) -> &dyn Any {
        self
    }
    #[inline(always)]
    fn as_any_mut(&mut self) -> &mut dyn Any {
        self
    }
    #[inline(always)]
    fn as_boxed_any(self: Box<Self>) -> Box<dyn Any> {
        self
    }
    #[inline(always)]
    fn type_name(&self) -> &'static str {
        type_name::<T>()
    }
    #[inline(always)]
    fn clone_object(&self) -> Box<dyn Variant> {
        Box::new(self.clone()) as Box<dyn Variant>
    }
}

impl dyn Variant {
    /// Is this [`Variant`] a specific type?
    #[inline(always)]
    #[must_use]
    pub fn is<T: Any>(&self) -> bool {
        TypeId::of::<T>() == self.type_id()
    }
}
