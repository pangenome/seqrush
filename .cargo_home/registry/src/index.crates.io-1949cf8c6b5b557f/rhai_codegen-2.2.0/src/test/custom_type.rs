#[cfg(not(feature = "metadata"))]
#[cfg(test)]
mod custom_type_tests {
    use crate::test::assert_streams_eq;
    use quote::quote;

    #[test]
    fn test_custom_type_tuple_struct() {
        let input = quote! {
            #[derive(Clone, CustomType)]
            pub struct Bar(
                #[rhai_type(skip)]
                #[cfg(not(feature = "no_float"))]
                rhai::FLOAT,
                INT,
                #[rhai_type(name = "boo", readonly)]
                String,
                Vec<INT>
            );
        };

        let result = crate::custom_type::derive_custom_type_impl(
            syn::parse2::<syn::DeriveInput>(input).unwrap(),
        );

        let expected = quote! {
            impl CustomType for Bar {
                fn build(mut builder: TypeBuilder<Self>) {
                    builder.with_name(stringify!(Bar));
                    builder.with_get_set("field1",
                        |obj: &mut Self| obj.1.clone(),
                        |obj: &mut Self, val| obj.1 = val
                    );
                    builder.with_get("boo", |obj: &mut Self| obj.2.clone());
                    builder.with_get_set("field3",
                        |obj: &mut Self| obj.3.clone(),
                        |obj: &mut Self, val| obj.3 = val
                    );
                }
            }
        };

        assert_streams_eq(result, expected);
    }

    #[test]
    fn test_custom_type_struct() {
        let input = quote! {
            #[derive(CustomType)]
            #[rhai_type(skip, name = "MyFoo", extra = Self::build_extra)]
            pub struct Foo {
                #[cfg(not(feature = "no_float"))]
                #[rhai_type(skip)]
                _dummy: rhai::FLOAT,
                #[rhai_type(get = get_bar)]
                pub bar: INT,
                #[rhai_type(name = "boo", readonly)]
                pub(crate) baz: String,
                #[rhai_type(set = Self::set_qux)]
                pub qux: Vec<INT>
            }
        };

        let result = crate::custom_type::derive_custom_type_impl(
            syn::parse2::<syn::DeriveInput>(input).unwrap(),
        );

        let expected = quote! {
            impl CustomType for Foo {
                fn build(mut builder: TypeBuilder<Self>) {
                    builder.with_name("MyFoo");
                    builder.with_get_set(stringify!(bar),
                        |obj: &mut Self| get_bar(&*obj),
                        |obj: &mut Self, val| obj.bar = val
                    );
                    builder.with_get("boo", |obj: &mut Self| obj.baz.clone());
                    builder.with_get_set(stringify!(qux),
                        |obj: &mut Self| obj.qux.clone(),
                        Self::set_qux
                    );
                    Self::build_extra(&mut builder);
                }
            }
        };

        assert_streams_eq(result, expected);
    }
}

#[cfg(feature = "metadata")]
#[cfg(test)]
mod custom_type_tests {
    use crate::test::assert_streams_eq;
    use quote::quote;

    #[test]
    fn test_custom_type_tuple_struct() {
        let input = quote! {
            /// Bar comments.
            #[derive(Clone, CustomType)]
            pub struct Bar(
                #[rhai_type(skip)]
                #[cfg(not(feature = "no_float"))]
                rhai::FLOAT,
                INT,
                /// boo comments.
                #[rhai_type(name = "boo", readonly)]
                String,
                /// This is a vector.
                Vec<INT>
            );
        };

        let result = crate::custom_type::derive_custom_type_impl(
            syn::parse2::<syn::DeriveInput>(input).unwrap(),
        );

        let expected = quote! {
            impl CustomType for Bar {
                fn build(mut builder: TypeBuilder<Self>) {
                    builder.with_name(stringify!(Bar)).with_comments(&"/// Bar comments.".lines().collect::<Vec<_>>()[..]);
                    builder.with_get_set("field1",
                        |obj: &mut Self| obj.1.clone(),
                        |obj: &mut Self, val| obj.1 = val
                    ).and_comments(&"".lines().collect::<Vec<_>>()[..]);
                    builder.with_get("boo", |obj: &mut Self| obj.2.clone())
                    .and_comments(&"/// boo comments.".lines().collect::<Vec<_>>()[..]);
                    builder.with_get_set("field3",
                        |obj: &mut Self| obj.3.clone(),
                        |obj: &mut Self, val| obj.3 = val
                    ).and_comments(&"/// This is a vector.".lines().collect::<Vec<_>>()[..]);
                }
            }
        };

        assert_streams_eq(result, expected);
    }

    #[test]
    fn test_custom_type_struct() {
        let input = quote! {
            /// Foo comments.
            #[derive(CustomType)]
            #[rhai_type(skip, name = "MyFoo", extra = Self::build_extra)]
            pub struct Foo {
                #[cfg(not(feature = "no_float"))]
                #[rhai_type(skip)]
                _dummy: rhai::FLOAT,
                #[rhai_type(get = get_bar)]
                pub bar: INT,
                /// boo comments.
                #[rhai_type(name = "boo", readonly)]
                pub(crate) baz: String,
                #[rhai_type(set = Self::set_qux)]
                pub qux: Vec<INT>
            }
        };

        let result = crate::custom_type::derive_custom_type_impl(
            syn::parse2::<syn::DeriveInput>(input).unwrap(),
        );

        let expected = quote! {
            impl CustomType for Foo {
                fn build(mut builder: TypeBuilder<Self>) {
                    builder.with_name("MyFoo").with_comments(&"/// Foo comments.".lines().collect::<Vec<_>>()[..]);
                    builder.with_get_set(stringify!(bar),
                        |obj: &mut Self| get_bar(&*obj),
                        |obj: &mut Self, val| obj.bar = val
                    ).and_comments(&"".lines().collect::<Vec<_>>()[..]);
                    builder.with_get("boo", |obj: &mut Self| obj.baz.clone())
                    .and_comments(&"/// boo comments.".lines().collect::<Vec<_>>()[..]);
                    builder.with_get_set(stringify!(qux),
                        |obj: &mut Self| obj.qux.clone(),
                        Self::set_qux
                    ).and_comments(&"".lines().collect::<Vec<_>>()[..]);
                    Self::build_extra(&mut builder);
                }
            }
        };

        assert_streams_eq(result, expected);
    }
}
