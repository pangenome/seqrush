//! Module implementing custom syntax for [`Engine`].
#![cfg(not(feature = "no_custom_syntax"))]

use crate::ast::Expr;
use crate::func::SendSync;
use crate::parser::ParseResult;
use crate::tokenizer::{is_reserved_keyword_or_symbol, is_valid_identifier, Token};
use crate::types::dynamic::Variant;
use crate::{
    Dynamic, Engine, EvalContext, Identifier, ImmutableString, LexError, Position, RhaiResult,
};
#[cfg(feature = "no_std")]
use std::prelude::v1::*;
use std::{borrow::Borrow, ops::Deref};

/// Collection of special markers for custom syntax definition.
pub mod markers {
    /// Special marker for matching an expression.
    pub const CUSTOM_SYNTAX_MARKER_EXPR: &str = "$expr$";
    /// Special marker for matching a statements block.
    pub const CUSTOM_SYNTAX_MARKER_BLOCK: &str = "$block$";
    /// Special marker for matching a function body.
    #[cfg(not(feature = "no_function"))]
    pub const CUSTOM_SYNTAX_MARKER_FUNC: &str = "$func$";
    /// Special marker for matching an identifier.
    pub const CUSTOM_SYNTAX_MARKER_IDENT: &str = "$ident$";
    /// Special marker for matching a single symbol.
    pub const CUSTOM_SYNTAX_MARKER_SYMBOL: &str = "$symbol$";
    /// Special marker for matching a string literal.
    pub const CUSTOM_SYNTAX_MARKER_STRING: &str = "$string$";
    /// Special marker for matching an integer number.
    pub const CUSTOM_SYNTAX_MARKER_INT: &str = "$int$";
    /// Special marker for matching a floating-point number.
    #[cfg(not(feature = "no_float"))]
    pub const CUSTOM_SYNTAX_MARKER_FLOAT: &str = "$float$";
    /// Special marker for matching a boolean value.
    pub const CUSTOM_SYNTAX_MARKER_BOOL: &str = "$bool$";
    /// Special marker for identifying the custom syntax variant.
    pub const CUSTOM_SYNTAX_MARKER_SYNTAX_VARIANT: &str = "$$";
}

/// A general expression evaluation trait object.
#[cfg(not(feature = "sync"))]
pub type FnCustomSyntaxEval = dyn Fn(&mut EvalContext, &[Expression], &Dynamic) -> RhaiResult;
/// A general expression evaluation trait object.
#[cfg(feature = "sync")]
pub type FnCustomSyntaxEval =
    dyn Fn(&mut EvalContext, &[Expression], &Dynamic) -> RhaiResult + Send + Sync;

/// A general expression parsing trait object.
#[cfg(not(feature = "sync"))]
pub type FnCustomSyntaxParse =
    dyn Fn(&[ImmutableString], &str, &mut Dynamic) -> ParseResult<Option<ImmutableString>>;
/// A general expression parsing trait object.
#[cfg(feature = "sync")]
pub type FnCustomSyntaxParse = dyn Fn(&[ImmutableString], &str, &mut Dynamic) -> ParseResult<Option<ImmutableString>>
    + Send
    + Sync;

/// An expression sub-tree in an [`AST`][crate::AST].
#[derive(Debug, Clone)]
pub struct Expression<'a>(&'a Expr);

impl<'a> From<&'a Expr> for Expression<'a> {
    #[inline(always)]
    fn from(expr: &'a Expr) -> Self {
        Self(expr)
    }
}

impl Expression<'_> {
    /// Evaluate this [expression tree][Expression] within an [evaluation context][`EvalContext`].
    ///
    /// # WARNING - Low Level API
    ///
    /// This function is very low level.  It evaluates an expression from an [`AST`][crate::AST].
    #[inline(always)]
    pub fn eval_with_context(&self, context: &mut EvalContext) -> RhaiResult {
        context.eval_expression_tree(self)
    }
    /// Evaluate this [expression tree][Expression] within an [evaluation context][`EvalContext`].
    ///
    /// The following option is available:
    ///
    /// * whether to rewind the [`Scope`][crate::Scope] after evaluation if the expression is a [`StmtBlock`][crate::ast::StmtBlock]
    ///
    /// # WARNING - Unstable API
    ///
    /// This API is volatile and may change in the future.
    ///
    /// # WARNING - Low Level API
    ///
    /// This function is _extremely_ low level.  It evaluates an expression from an [`AST`][crate::AST].
    #[deprecated = "This API is NOT deprecated, but it is considered volatile and may change in the future."]
    #[inline(always)]
    pub fn eval_with_context_raw(
        &self,
        context: &mut EvalContext,
        rewind_scope: bool,
    ) -> RhaiResult {
        #[allow(deprecated)]
        context.eval_expression_tree_raw(self, rewind_scope)
    }
    /// Get the value of this expression if it is a variable name or a string constant.
    ///
    /// Returns [`None`] also if the constant is not of the specified type.
    #[inline(always)]
    #[must_use]
    pub fn get_string_value(&self) -> Option<&str> {
        match self.0 {
            #[cfg(not(feature = "no_module"))]
            Expr::Variable(x, ..) if !x.2.is_empty() => None,
            Expr::Variable(x, ..) => Some(&x.1),
            #[cfg(not(feature = "no_function"))]
            Expr::ThisPtr(..) => Some(crate::engine::KEYWORD_THIS),
            Expr::StringConstant(x, ..) => Some(x),
            _ => None,
        }
    }
    /// Get the position of this expression.
    #[inline(always)]
    #[must_use]
    pub const fn position(&self) -> Position {
        self.0.position()
    }
    /// Get the value of this expression if it is a literal constant.
    ///
    /// Supports [`INT`][crate::INT], [`FLOAT`][crate::FLOAT], `()`, `char`, `bool` and
    /// [`ImmutableString`][crate::ImmutableString].
    ///
    /// Returns [`None`] also if the constant is not of the specified type.
    #[inline]
    #[must_use]
    pub fn get_literal_value<T: Variant + Clone>(&self) -> Option<T> {
        // Coded this way in order to maximally leverage potentials for dead-code removal.
        match self.0 {
            Expr::DynamicConstant(x, ..) => x.clone().try_cast::<T>(),
            Expr::IntegerConstant(x, ..) => reify! { *x => Option<T> },

            #[cfg(not(feature = "no_float"))]
            Expr::FloatConstant(x, ..) => reify! { *x => Option<T> },

            Expr::CharConstant(x, ..) => reify! { *x => Option<T> },
            Expr::StringConstant(x, ..) => reify! { x.clone() => Option<T> },
            Expr::Variable(x, ..) => reify! { x.1.clone() => Option<T> },
            Expr::BoolConstant(x, ..) => reify! { *x => Option<T> },
            Expr::Unit(..) => reify! { () => Option<T> },

            _ => None,
        }
    }
}

impl Borrow<Expr> for Expression<'_> {
    #[inline(always)]
    fn borrow(&self) -> &Expr {
        self.0
    }
}

impl AsRef<Expr> for Expression<'_> {
    #[inline(always)]
    fn as_ref(&self) -> &Expr {
        self.0
    }
}

impl Deref for Expression<'_> {
    type Target = Expr;

    #[inline(always)]
    fn deref(&self) -> &Self::Target {
        self.0
    }
}

/// Definition of a custom syntax definition.
pub struct CustomSyntax {
    /// A parsing function to return the next token in a custom syntax based on the
    /// symbols parsed so far.
    pub parse: Box<FnCustomSyntaxParse>,
    /// Custom syntax implementation function.
    pub func: Box<FnCustomSyntaxEval>,
    /// Any variables added/removed in the scope?
    pub scope_may_be_changed: bool,
}

impl Engine {
    /// Register a custom syntax with the [`Engine`].
    ///
    /// Not available under `no_custom_syntax`.
    ///
    /// * `symbols` holds a slice of strings that define the custom syntax.
    /// * `scope_may_be_changed` specifies variables _may_ be added/removed by this custom syntax.
    /// * `func` is the implementation function.
    ///
    /// ## Note on `symbols`
    ///
    /// * Whitespaces around symbols are stripped.
    /// * Symbols that are all-whitespace or empty are ignored.
    /// * If `symbols` does not contain at least one valid token, then the custom syntax registration
    ///   is simply ignored.
    ///
    /// ## Note on `scope_may_be_changed`
    ///
    /// If `scope_may_be_changed` is `true`, then _size_ of the current [`Scope`][crate::Scope]
    /// _may_ be modified by this custom syntax.
    ///
    /// Adding new variables and/or removing variables count.
    ///
    /// Simply modifying the values of existing variables does NOT count, as the _size_ of the
    /// current [`Scope`][crate::Scope] is unchanged, so `false` should be passed.
    ///
    /// Replacing one variable with another (i.e. adding a new variable and removing one variable at
    /// the same time so that the total _size_ of the [`Scope`][crate::Scope] is unchanged) also
    /// does NOT count, so `false` should be passed.
    pub fn register_custom_syntax<S: AsRef<str> + Into<Identifier>>(
        &mut self,
        symbols: impl AsRef<[S]>,
        scope_may_be_changed: bool,
        func: impl Fn(&mut EvalContext, &[Expression]) -> RhaiResult + SendSync + 'static,
    ) -> ParseResult<&mut Self> {
        #[allow(clippy::wildcard_imports)]
        use markers::*;

        let mut segments = Vec::<ImmutableString>::new();

        for s in symbols.as_ref() {
            let s = s.as_ref().trim();

            // Skip empty symbols
            if s.is_empty() {
                continue;
            }

            let token = Token::lookup_symbol_from_syntax(s).or_else(|| {
                is_reserved_keyword_or_symbol(s)
                    .0
                    .then(|| Token::Reserved(Box::new(s.into())))
            });

            let seg = match s {
                // Markers not in first position
                CUSTOM_SYNTAX_MARKER_IDENT
                | CUSTOM_SYNTAX_MARKER_SYMBOL
                | CUSTOM_SYNTAX_MARKER_EXPR
                | CUSTOM_SYNTAX_MARKER_BLOCK
                | CUSTOM_SYNTAX_MARKER_BOOL
                | CUSTOM_SYNTAX_MARKER_INT
                | CUSTOM_SYNTAX_MARKER_STRING
                    if !segments.is_empty() =>
                {
                    s.into()
                }

                // Markers not in first position
                #[cfg(not(feature = "no_function"))]
                CUSTOM_SYNTAX_MARKER_FUNC if !segments.is_empty() => s.into(),

                // Markers not in first position
                #[cfg(not(feature = "no_float"))]
                CUSTOM_SYNTAX_MARKER_FLOAT if !segments.is_empty() => s.into(),

                // Identifier not in first position
                _ if !segments.is_empty() && is_valid_identifier(s) => s.into(),

                // Keyword/symbol not in first position
                _ if !segments.is_empty() && token.is_some() => {
                    // Make it a custom keyword/symbol if it is disabled or reserved
                    if (self.is_symbol_disabled(s)
                        || token.as_ref().map_or(false, Token::is_reserved))
                        && !self.custom_keywords.contains_key(s)
                    {
                        self.custom_keywords.insert(s.into(), None);
                    }
                    s.into()
                }

                // Standard keyword in first position but not disabled
                _ if segments.is_empty()
                    && token.as_ref().map_or(false, Token::is_standard_keyword)
                    && !self.is_symbol_disabled(s) =>
                {
                    return Err(LexError::ImproperSymbol(
                        s.to_string(),
                        format!("Improper symbol for custom syntax at position #0: '{s}'"),
                    )
                    .into_err(Position::NONE));
                }

                // Identifier or symbol in first position
                _ if segments.is_empty()
                    && (is_valid_identifier(s) || is_reserved_keyword_or_symbol(s).0) =>
                {
                    // Make it a custom keyword/symbol if it is disabled or reserved
                    if self.is_symbol_disabled(s)
                        || (token.as_ref().map_or(false, Token::is_reserved)
                            && !self.custom_keywords.contains_key(s))
                    {
                        self.custom_keywords.insert(s.into(), None);
                    }
                    s.into()
                }

                // Anything else is an error
                _ => {
                    return Err(LexError::ImproperSymbol(
                        s.to_string(),
                        format!(
                            "Improper symbol for custom syntax at position #{}: '{s}'",
                            segments.len() + 1,
                        ),
                    )
                    .into_err(Position::NONE));
                }
            };

            segments.push(seg);
        }

        // If the syntax has nothing, just ignore the registration
        if segments.is_empty() {
            return Ok(self);
        }

        // The first keyword/symbol is the discriminator
        let key = segments[0].clone();

        self.register_custom_syntax_with_state_raw(
            key,
            // Construct the parsing function
            move |stream, _, _| match stream.len() {
                len if len >= segments.len() => Ok(None),
                len => Ok(Some(segments[len].clone())),
            },
            scope_may_be_changed,
            move |context, expressions, _| func(context, expressions),
        );

        Ok(self)
    }
    /// Register a custom syntax with the [`Engine`] with custom user-defined state.
    ///
    /// Not available under `no_custom_syntax`.
    ///
    /// # WARNING - Low Level API
    ///
    /// This function is very low level.
    ///
    /// * `scope_may_be_changed` specifies variables have been added/removed by this custom syntax.
    /// * `parse` is the parsing function.
    /// * `func` is the implementation function.
    ///
    /// All custom keywords used as symbols must be manually registered via [`Engine::register_custom_operator`].
    /// Otherwise, they won't be recognized.
    ///
    /// # Parsing Function Signature
    ///
    /// The parsing function has the following signature:
    ///
    /// `Fn(symbols: &[ImmutableString], look_ahead: &str, state: &mut Dynamic) -> Result<Option<ImmutableString>, ParseError>`
    ///
    /// where:
    /// * `symbols`: a slice of symbols that have been parsed so far, possibly containing `$expr$` and/or `$block$`;
    ///   `$ident$` and other literal markers are replaced by the actual text
    /// * `look_ahead`: a string slice containing the next symbol that is about to be read
    /// * `state`: a [`Dynamic`] value that contains a user-defined state
    ///
    /// ## Return value
    ///
    /// * `Ok(None)`: parsing complete and there are no more symbols to match.
    /// * `Ok(Some(symbol))`: the next symbol to match, which can also be `$expr$`, `$ident$` or `$block$`.
    /// * `Err(ParseError)`: error that is reflected back to the [`Engine`], normally `ParseError(ParseErrorType::BadInput(LexError::ImproperSymbol(message)), Position::NONE)` to indicate a syntax error, but it can be any [`ParseError`][crate::ParseError].
    pub fn register_custom_syntax_with_state_raw(
        &mut self,
        key: impl Into<Identifier>,
        parse: impl Fn(&[ImmutableString], &str, &mut Dynamic) -> ParseResult<Option<ImmutableString>>
            + SendSync
            + 'static,
        scope_may_be_changed: bool,
        func: impl Fn(&mut EvalContext, &[Expression], &Dynamic) -> RhaiResult + SendSync + 'static,
    ) -> &mut Self {
        self.custom_syntax.insert(
            key.into(),
            CustomSyntax {
                parse: Box::new(parse),
                func: Box::new(func),
                scope_may_be_changed,
            }
            .into(),
        );
        self
    }
}
