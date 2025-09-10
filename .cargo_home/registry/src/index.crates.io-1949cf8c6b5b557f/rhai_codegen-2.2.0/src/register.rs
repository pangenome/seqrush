use quote::{quote, quote_spanned};
use syn::{parse::Parser, spanned::Spanned};

pub fn generated_module_path(
    fn_path: &syn::Path,
) -> syn::punctuated::Punctuated<syn::PathSegment, syn::Token![::]> {
    let mut g = fn_path.clone().segments;
    g.pop();
    let ident = syn::Ident::new(
        &format!("rhai_fn_{}", fn_path.segments.last().unwrap().ident),
        fn_path.span(),
    );
    g.push_value(syn::PathSegment {
        ident,
        arguments: syn::PathArguments::None,
    });
    g
}

type RegisterMacroInput = (syn::Expr, proc_macro2::TokenStream, syn::Path);

pub fn parse_register_macro(
    args: proc_macro::TokenStream,
) -> Result<RegisterMacroInput, syn::Error> {
    let args = syn::punctuated::Punctuated::<_, syn::Token![,]>::parse_separated_nonempty
        .parse(args)
        .unwrap();
    let arg_span = args.span();
    let mut items: Vec<_> = args.into_iter().collect();
    if items.len() != 3 {
        return Err(syn::Error::new(
            arg_span,
            "this macro requires three arguments",
        ));
    }
    let export_name = match &items[1] {
        syn::Expr::Lit(lit_str) => quote_spanned!(items[1].span() => #lit_str),
        expr => quote! { #expr },
    };
    let rust_mod_path = match &items[2] {
        syn::Expr::Path(ref path) => path.path.clone(),
        _ => {
            return Err(syn::Error::new(
                items[2].span(),
                "third argument must be a function name",
            ))
        }
    };
    let module = items.remove(0);
    Ok((module, export_name, rust_mod_path))
}
