#![cfg(not(feature = "no_std"))]
#![cfg(any(not(target_family = "wasm"), not(target_os = "unknown")))]

use crate::eval::GlobalRuntimeState;
use crate::func::{locked_read, locked_write};
use crate::{
    Engine, Identifier, Locked, Module, ModuleResolver, Position, RhaiResultOf, Scope, Shared,
    SharedModule, ERR,
};

use std::{
    collections::BTreeMap,
    io::Error as IoError,
    path::{Path, PathBuf},
};

pub const RHAI_SCRIPT_EXTENSION: &str = "rhai";

/// A [module][Module] resolution service that loads [module][Module] script files from the file system.
///
/// ## Caching
///
/// Resolved [Modules][Module] are cached internally so script files are not reloaded and recompiled
/// for subsequent requests.
///
/// Use [`clear_cache`][FileModuleResolver::clear_cache] or
/// [`clear_cache_for_path`][FileModuleResolver::clear_cache_for_path] to clear the internal cache.
///
/// ## Namespace
///
/// When a function within a script file module is called, all functions defined within the same
/// script are available, evan `private` ones.  In other words, functions defined in a module script
/// can always cross-call each other.
///
/// # Example
///
/// ```
/// use rhai::Engine;
/// use rhai::module_resolvers::FileModuleResolver;
///
/// // Create a new 'FileModuleResolver' loading scripts from the 'scripts' subdirectory
/// // with file extension '.x'.
/// let resolver = FileModuleResolver::new_with_path_and_extension("./scripts", "x");
///
/// let mut engine = Engine::new();
///
/// engine.set_module_resolver(resolver);
/// ```
#[derive(Debug)]
pub struct FileModuleResolver {
    /// Base path of the directory holding script files.
    base_path: Option<PathBuf>,
    /// File extension of script files, default `.rhai`.
    extension: Identifier,
    /// Is the cache enabled?
    cache_enabled: bool,
    /// [`Scope`] holding variables for compiling scripts.
    scope: Scope<'static>,
    /// Internal cache of resolved modules.
    ///
    /// The cache is wrapped in interior mutability because [`resolve`][FileModuleResolver::resolve]
    /// is immutable.
    cache: Locked<BTreeMap<PathBuf, SharedModule>>,
}

impl Default for FileModuleResolver {
    #[inline(always)]
    fn default() -> Self {
        Self::new()
    }
}

impl FileModuleResolver {
    /// Create a new [`FileModuleResolver`] with the current directory as base path.
    ///
    /// The default extension is `.rhai`.
    ///
    /// # Example
    ///
    /// ```
    /// use rhai::Engine;
    /// use rhai::module_resolvers::FileModuleResolver;
    ///
    /// // Create a new 'FileModuleResolver' loading scripts from the current directory
    /// // with file extension '.rhai' (the default).
    /// let resolver = FileModuleResolver::new();
    ///
    /// let mut engine = Engine::new();
    /// engine.set_module_resolver(resolver);
    /// ```
    #[inline(always)]
    #[must_use]
    pub fn new() -> Self {
        Self::new_with_extension(RHAI_SCRIPT_EXTENSION)
    }

    /// Create a new [`FileModuleResolver`] with a specific base path.
    ///
    /// The default extension is `.rhai`.
    ///
    /// # Example
    ///
    /// ```
    /// use rhai::Engine;
    /// use rhai::module_resolvers::FileModuleResolver;
    ///
    /// // Create a new 'FileModuleResolver' loading scripts from the 'scripts' subdirectory
    /// // with file extension '.rhai' (the default).
    /// let resolver = FileModuleResolver::new_with_path("./scripts");
    ///
    /// let mut engine = Engine::new();
    /// engine.set_module_resolver(resolver);
    /// ```
    #[inline(always)]
    #[must_use]
    pub fn new_with_path(path: impl Into<PathBuf>) -> Self {
        Self::new_with_path_and_extension(path, RHAI_SCRIPT_EXTENSION)
    }

    /// Create a new [`FileModuleResolver`] with a file extension.
    ///
    /// # Example
    ///
    /// ```
    /// use rhai::Engine;
    /// use rhai::module_resolvers::FileModuleResolver;
    ///
    /// // Create a new 'FileModuleResolver' loading scripts with file extension '.rhai' (the default).
    /// let resolver = FileModuleResolver::new_with_extension("rhai");
    ///
    /// let mut engine = Engine::new();
    /// engine.set_module_resolver(resolver);
    /// ```
    #[inline(always)]
    #[must_use]
    pub fn new_with_extension(extension: impl Into<Identifier>) -> Self {
        Self {
            base_path: None,
            extension: extension.into(),
            cache_enabled: true,
            cache: BTreeMap::new().into(),
            scope: Scope::new(),
        }
    }

    /// Create a new [`FileModuleResolver`] with a specific base path and file extension.
    ///
    /// # Example
    ///
    /// ```
    /// use rhai::Engine;
    /// use rhai::module_resolvers::FileModuleResolver;
    ///
    /// // Create a new 'FileModuleResolver' loading scripts from the 'scripts' subdirectory
    /// // with file extension '.x'.
    /// let resolver = FileModuleResolver::new_with_path_and_extension("./scripts", "x");
    ///
    /// let mut engine = Engine::new();
    /// engine.set_module_resolver(resolver);
    /// ```
    #[inline(always)]
    #[must_use]
    pub fn new_with_path_and_extension(
        path: impl Into<PathBuf>,
        extension: impl Into<Identifier>,
    ) -> Self {
        Self {
            base_path: Some(path.into()),
            extension: extension.into(),
            cache_enabled: true,
            cache: BTreeMap::new().into(),
            scope: Scope::new(),
        }
    }

    /// Get the base path for script files.
    #[inline(always)]
    #[must_use]
    pub fn base_path(&self) -> Option<&Path> {
        self.base_path.as_deref()
    }
    /// Set the base path for script files.
    #[inline(always)]
    pub fn set_base_path(&mut self, path: impl Into<PathBuf>) -> &mut Self {
        self.base_path = Some(path.into());
        self
    }

    /// Get the script file extension.
    #[inline(always)]
    #[must_use]
    pub fn extension(&self) -> &str {
        &self.extension
    }

    /// Set the script file extension.
    #[inline(always)]
    pub fn set_extension(&mut self, extension: impl Into<Identifier>) -> &mut Self {
        self.extension = extension.into();
        self
    }

    /// Get a reference to the file module resolver's [scope][Scope].
    ///
    /// The [scope][Scope] is used for compiling module scripts.
    #[inline(always)]
    #[must_use]
    pub const fn scope(&self) -> &Scope {
        &self.scope
    }

    /// Set the file module resolver's [scope][Scope].
    ///
    /// The [scope][Scope] is used for compiling module scripts.
    #[inline(always)]
    pub fn set_scope(&mut self, scope: Scope<'static>) {
        self.scope = scope;
    }

    /// Get a mutable reference to the file module resolver's [scope][Scope].
    ///
    /// The [scope][Scope] is used for compiling module scripts.
    #[inline(always)]
    #[must_use]
    pub fn scope_mut(&mut self) -> &mut Scope<'static> {
        &mut self.scope
    }

    /// Enable/disable the cache.
    #[inline(always)]
    pub fn enable_cache(&mut self, enable: bool) -> &mut Self {
        self.cache_enabled = enable;
        self
    }
    /// Is the cache enabled?
    #[inline(always)]
    #[must_use]
    pub const fn is_cache_enabled(&self) -> bool {
        self.cache_enabled
    }

    /// Is a particular path cached?
    #[inline]
    #[must_use]
    pub fn is_cached(&self, path: impl AsRef<Path>) -> bool {
        if !self.cache_enabled {
            return false;
        }
        locked_read(&self.cache)
            .unwrap()
            .contains_key(path.as_ref())
    }
    /// Empty the internal cache.
    #[inline]
    pub fn clear_cache(&mut self) -> &mut Self {
        locked_write(&self.cache).unwrap().clear();
        self
    }
    /// Remove the specified path from internal cache.
    ///
    /// The next time this path is resolved, the script file will be loaded once again.
    #[inline]
    #[must_use]
    pub fn clear_cache_for_path(&mut self, path: impl AsRef<Path>) -> Option<SharedModule> {
        locked_write(&self.cache)
            .unwrap()
            .remove_entry(path.as_ref())
            .map(|(.., v)| v)
    }
    /// Construct a full file path.
    #[must_use]
    pub fn get_file_path(&self, path: &str, source_path: Option<&Path>) -> PathBuf {
        let path = Path::new(path);

        let mut file_path;

        if path.is_relative() {
            file_path = self
                .base_path
                .clone()
                .or_else(|| source_path.map(Into::into))
                .unwrap_or_default();
            file_path.push(path);
        } else {
            file_path = path.into();
        }

        file_path.set_extension(self.extension.as_str()); // Force extension
        file_path
    }

    /// Resolve a module based on a path.
    fn impl_resolve(
        &self,
        engine: &Engine,
        global: &mut GlobalRuntimeState,
        scope: &mut Scope,
        source: Option<&str>,
        path: &str,
        pos: Position,
    ) -> Result<SharedModule, Box<crate::EvalAltResult>> {
        // Load relative paths from source if there is no base path specified
        let source_path = global
            .source()
            .or(source)
            .and_then(|p| Path::new(p).parent());

        let file_path = self.get_file_path(path, source_path);

        if self.is_cache_enabled() {
            if let Some(module) = locked_read(&self.cache).unwrap().get(&file_path) {
                return Ok(module.clone());
            }
        }

        let mut ast = engine
            .compile_file_with_scope(&self.scope, file_path.clone())
            .map_err(|err| match *err {
                ERR::ErrorSystem(.., err) if err.is::<IoError>() => {
                    Box::new(ERR::ErrorModuleNotFound(path.to_string(), pos))
                }
                _ => Box::new(ERR::ErrorInModule(path.to_string(), err, pos)),
            })?;

        ast.set_source(path);

        let m: Shared<_> = Module::eval_ast_as_new_raw(engine, scope, global, &ast)
            .map_err(|err| Box::new(ERR::ErrorInModule(path.to_string(), err, pos)))?
            .into();

        if self.is_cache_enabled() {
            locked_write(&self.cache)
                .unwrap()
                .insert(file_path, m.clone());
        }

        Ok(m)
    }
}

impl ModuleResolver for FileModuleResolver {
    fn resolve_raw(
        &self,
        engine: &Engine,
        global: &mut GlobalRuntimeState,
        scope: &mut Scope,
        path: &str,
        pos: Position,
    ) -> RhaiResultOf<SharedModule> {
        self.impl_resolve(engine, global, scope, None, path, pos)
    }

    #[inline(always)]
    fn resolve(
        &self,
        engine: &Engine,
        source: Option<&str>,
        path: &str,
        pos: Position,
    ) -> RhaiResultOf<SharedModule> {
        let global = &mut engine.new_global_runtime_state();
        let scope = &mut Scope::new();
        self.impl_resolve(engine, global, scope, source, path, pos)
    }

    /// Resolve an `AST` based on a path string.
    ///
    /// The file system is accessed during each call; the internal cache is by-passed.
    fn resolve_ast(
        &self,
        engine: &Engine,
        source_path: Option<&str>,
        path: &str,
        pos: Position,
    ) -> Option<RhaiResultOf<crate::AST>> {
        // Construct the script file path
        let file_path = self.get_file_path(path, source_path.map(Path::new));

        // Load the script file and compile it
        Some(
            engine
                .compile_file(file_path)
                .map(|mut ast| {
                    ast.set_source(path);
                    ast
                })
                .map_err(|err| match *err {
                    ERR::ErrorSystem(.., err) if err.is::<IoError>() => {
                        ERR::ErrorModuleNotFound(path.to_string(), pos).into()
                    }
                    _ => ERR::ErrorInModule(path.to_string(), err, pos).into(),
                }),
        )
    }
}
