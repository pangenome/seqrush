use crate::{
    Engine, Identifier, Module, ModuleResolver, Position, RhaiResultOf, SharedModule, SmartString,
    ERR,
};
#[cfg(feature = "no_std")]
use std::prelude::v1::*;
use std::{
    collections::btree_map::{IntoIter, Iter},
    collections::BTreeMap,
    ops::AddAssign,
};

/// A static [module][Module] resolution service that serves [modules][Module] added into it.
///
/// # Example
///
/// ```
/// use rhai::{Engine, Module};
/// use rhai::module_resolvers::StaticModuleResolver;
///
/// let mut resolver = StaticModuleResolver::new();
///
/// let module = Module::new();
/// resolver.insert("hello", module);
///
/// let mut engine = Engine::new();
///
/// engine.set_module_resolver(resolver);
/// ```
#[derive(Debug, Clone, Default)]
pub struct StaticModuleResolver(BTreeMap<Identifier, SharedModule>);

impl StaticModuleResolver {
    /// Create a new [`StaticModuleResolver`].
    ///
    /// # Example
    ///
    /// ```
    /// use rhai::{Engine, Module};
    /// use rhai::module_resolvers::StaticModuleResolver;
    ///
    /// let mut resolver = StaticModuleResolver::new();
    ///
    /// let module = Module::new();
    /// resolver.insert("hello", module);
    ///
    /// let mut engine = Engine::new();
    /// engine.set_module_resolver(resolver);
    /// ```
    #[inline(always)]
    #[must_use]
    pub const fn new() -> Self {
        Self(BTreeMap::new())
    }
    /// Add a [module][Module] keyed by its path.
    #[inline]
    pub fn insert(&mut self, path: impl Into<Identifier>, mut module: Module) {
        let path = path.into();

        if module.id().is_none() {
            module.set_id(path.clone());
        }

        module.build_index();
        self.0.insert(path, module.into());
    }
    /// Remove a [module][Module] given its path.
    #[inline(always)]
    pub fn remove(&mut self, path: &str) -> Option<SharedModule> {
        self.0.remove(path)
    }
    /// Does the path exist?
    #[inline(always)]
    #[must_use]
    pub fn contains_path(&self, path: &str) -> bool {
        self.0.contains_key(path)
    }
    /// Get an iterator of all the [modules][Module].
    #[inline]
    pub fn iter(&self) -> impl Iterator<Item = (&str, &SharedModule)> {
        self.0.iter().map(|(k, v)| (k.as_str(), v))
    }
    /// Get a mutable iterator of all the [modules][Module].
    #[inline]
    pub fn iter_mut(&mut self) -> impl Iterator<Item = (&str, &mut SharedModule)> {
        self.0.iter_mut().map(|(k, v)| (k.as_str(), v))
    }
    /// Get an iterator of all the [module][Module] paths.
    #[inline]
    pub fn paths(&self) -> impl Iterator<Item = &str> {
        self.0.keys().map(SmartString::as_str)
    }
    /// Get an iterator of all the [modules][Module].
    #[inline(always)]
    pub fn values(&self) -> impl Iterator<Item = &SharedModule> {
        self.0.values()
    }
    /// Remove all [modules][Module].
    #[inline(always)]
    pub fn clear(&mut self) -> &mut Self {
        self.0.clear();
        self
    }
    /// Returns `true` if this [`StaticModuleResolver`] contains no module resolvers.
    #[inline(always)]
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }
    /// Get the number of [modules][Module] in this [`StaticModuleResolver`].
    #[inline(always)]
    #[must_use]
    pub fn len(&self) -> usize {
        self.0.len()
    }
    /// Merge another [`StaticModuleResolver`] into this.
    /// The other [`StaticModuleResolver`] is consumed.
    ///
    /// Existing modules of the same path name are overwritten.
    #[inline]
    pub fn merge(&mut self, other: Self) -> &mut Self {
        self.0.extend(other.0);
        self
    }
}

impl IntoIterator for StaticModuleResolver {
    type Item = (Identifier, SharedModule);
    type IntoIter = IntoIter<SmartString, SharedModule>;

    #[inline(always)]
    fn into_iter(self) -> Self::IntoIter {
        self.0.into_iter()
    }
}

impl<'a> IntoIterator for &'a StaticModuleResolver {
    type Item = (&'a Identifier, &'a SharedModule);
    type IntoIter = Iter<'a, SmartString, SharedModule>;

    #[inline(always)]
    fn into_iter(self) -> Self::IntoIter {
        self.0.iter()
    }
}

impl AddAssign<Self> for StaticModuleResolver {
    #[inline(always)]
    fn add_assign(&mut self, rhs: Self) {
        self.merge(rhs);
    }
}

impl ModuleResolver for StaticModuleResolver {
    #[inline]
    fn resolve(
        &self,
        _: &Engine,
        _: Option<&str>,
        path: &str,
        pos: Position,
    ) -> RhaiResultOf<SharedModule> {
        self.0
            .get(path)
            .cloned()
            .ok_or_else(|| ERR::ErrorModuleNotFound(path.into(), pos).into())
    }
}
