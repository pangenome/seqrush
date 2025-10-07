/// Macro to cast an identifier or expression to another type with type checks.
///
/// Runs _code_ if _variable_ or _expression_ is of type _type_, otherwise run _fallback_.
///
/// # Syntax
///
/// * `reify! { `_variable_ or _expression_` => |`_temp-variable_`: `_type_`|` _code_`,` `||` _fallback_ `)`
/// * `reify! { `_variable_ or _expression_` => |`_temp-variable_`: `_type_`|` _code_ `)`
/// * `reify! { `_variable_ or _expression_ `=>` `Option<`_type_`>` `)`
/// * `reify! { `_variable_ or _expression_ `=>` _type_ `)`
///
/// * `reify! { `_expression_ `=> !!!` _type_ `)`  (unsafe, no type checks!)
macro_rules! reify {
    ($old:ident => |$new:ident : $t:ty| $code:expr, || $fallback:expr) => {{
        #[allow(clippy::redundant_else)]
        if std::any::TypeId::of::<$t>() == std::any::Any::type_id(&$old) {
            // SAFETY: This is safe because we already checked to make sure the two types
            // are actually the same.
            let $new: $t = unsafe { std::mem::transmute_copy(&std::mem::ManuallyDrop::new($old)) };
            $code
        } else {
            $fallback
        }
    }};
    ($old:expr => |$new:ident : $t:ty| $code:expr, || $fallback:expr) => {{
        let old = $old;
        reify! { old => |$new: $t| $code, || $fallback }
    }};

    ($old:ident => |$new:ident : $t:ty| $code:expr) => {
        reify! { $old => |$new: $t| $code, || () }
    };
    ($old:expr => |$new:ident : $t:ty| $code:expr) => {
        reify! { $old => |$new: $t| $code, || () }
    };

    ($old:ident => Option<$t:ty>) => {
        reify! { $old => |v: $t| Some(v), || None }
    };
    ($old:expr => Option<$t:ty>) => {
        reify! { $old => |v: $t| Some(v), || None }
    };

    ($old:ident => $t:ty) => {
        reify! { $old => |v: $t| v, || unreachable!() }
    };
    ($old:expr => $t:ty) => {
        reify! { $old => |v: $t| v, || unreachable!() }
    };

    ($old:expr => !!! $t:ty) => {{
        let old_value = $old;
        let new_value: $t =
            // SAFETY: This is really dangerous (thus the exclamation marks).
            //         Make sure that the transmute is safe.
            unsafe { std::mem::transmute_copy(&std::mem::ManuallyDrop::new(old_value)) };
        new_value
    }};
}
