//! Implementations of [`serde::Deserialize`].

use crate::{Dynamic, Identifier, ImmutableString, Scope, INT};
use serde::{
    de::{Error, SeqAccess, Visitor},
    Deserialize, Deserializer,
};
use std::fmt;
#[cfg(feature = "no_std")]
use std::prelude::v1::*;

#[cfg(feature = "decimal")]
use num_traits::FromPrimitive;

struct DynamicVisitor;

impl<'de> Visitor<'de> for DynamicVisitor {
    type Value = Dynamic;

    #[cold]
    #[inline(never)]
    fn expecting(&self, f: &mut fmt::Formatter) -> fmt::Result {
        f.write_str("any type that can be converted into a Dynamic")
    }
    #[inline(always)]
    fn visit_bool<E: Error>(self, v: bool) -> Result<Self::Value, E> {
        Ok(v.into())
    }
    #[inline(always)]
    fn visit_i8<E: Error>(self, v: i8) -> Result<Self::Value, E> {
        Ok(INT::from(v).into())
    }
    #[inline(always)]
    fn visit_i16<E: Error>(self, v: i16) -> Result<Self::Value, E> {
        Ok(INT::from(v).into())
    }
    #[inline(always)]
    fn visit_i32<E: Error>(self, v: i32) -> Result<Self::Value, E> {
        Ok(INT::from(v).into())
    }
    #[inline]
    fn visit_i64<E: Error>(self, v: i64) -> Result<Self::Value, E> {
        #[cfg(not(feature = "only_i32"))]
        return Ok(v.into());

        #[cfg(feature = "only_i32")]
        if v <= INT::MAX as i64 {
            return Ok(Dynamic::from(v as INT));
        }

        #[allow(unreachable_code)]
        {
            #[cfg(feature = "decimal")]
            if let Some(n) = rust_decimal::Decimal::from_i64(v) {
                return Ok(Dynamic::from_decimal(n));
            }

            #[cfg(not(feature = "no_float"))]
            return Ok(Dynamic::from_float(v as crate::FLOAT));

            Err(Error::custom(format!("integer number too large: {v}")))
        }
    }
    #[inline]
    fn visit_i128<E: Error>(self, v: i128) -> Result<Self::Value, E> {
        if v <= i128::from(INT::MAX) {
            return Ok(Dynamic::from(v as INT));
        }

        #[allow(unreachable_code)]
        {
            #[cfg(feature = "decimal")]
            if let Some(n) = rust_decimal::Decimal::from_i128(v) {
                return Ok(Dynamic::from_decimal(n));
            }

            #[cfg(not(feature = "no_float"))]
            return Ok(Dynamic::from_float(v as crate::FLOAT));

            Err(Error::custom(format!("integer number too large: {v}")))
        }
    }
    #[inline(always)]
    fn visit_u8<E: Error>(self, v: u8) -> Result<Self::Value, E> {
        Ok(INT::from(v).into())
    }
    #[inline(always)]
    fn visit_u16<E: Error>(self, v: u16) -> Result<Self::Value, E> {
        Ok(INT::from(v).into())
    }
    #[inline]
    fn visit_u32<E: Error>(self, v: u32) -> Result<Self::Value, E> {
        #[cfg(not(feature = "only_i32"))]
        return Ok(Dynamic::from(v as INT));

        #[cfg(feature = "only_i32")]
        if v <= INT::MAX as u32 {
            return Ok(Dynamic::from(v as INT));
        }

        #[allow(unreachable_code)]
        {
            #[cfg(feature = "decimal")]
            if let Some(n) = rust_decimal::Decimal::from_u32(v) {
                return Ok(Dynamic::from_decimal(n));
            }

            #[cfg(not(feature = "no_float"))]
            return Ok(Dynamic::from_float(v as crate::FLOAT));

            Err(Error::custom(format!("integer number too large: {v}")))
        }
    }
    #[inline]
    fn visit_u64<E: Error>(self, v: u64) -> Result<Self::Value, E> {
        if v <= INT::MAX as u64 {
            return Ok(Dynamic::from(v as INT));
        }

        #[cfg(feature = "decimal")]
        if let Some(n) = rust_decimal::Decimal::from_u64(v) {
            return Ok(Dynamic::from_decimal(n));
        }

        #[cfg(not(feature = "no_float"))]
        return Ok(Dynamic::from_float(v as crate::FLOAT));

        #[allow(unreachable_code)]
        Err(Error::custom(format!("integer number too large: {v}")))
    }
    #[inline]
    fn visit_u128<E: Error>(self, v: u128) -> Result<Self::Value, E> {
        if v <= INT::MAX as u128 {
            return Ok(Dynamic::from(v as INT));
        }

        #[cfg(feature = "decimal")]
        if let Some(n) = rust_decimal::Decimal::from_u128(v) {
            return Ok(Dynamic::from_decimal(n));
        }

        #[cfg(not(feature = "no_float"))]
        return Ok(Dynamic::from_float(v as crate::FLOAT));

        #[allow(unreachable_code)]
        Err(Error::custom(format!("integer number too large: {v}")))
    }

    #[cfg(not(feature = "no_float"))]
    #[inline(always)]
    fn visit_f32<E: Error>(self, v: f32) -> Result<Self::Value, E> {
        #[cfg(not(feature = "no_float"))]
        return Ok((v as crate::FLOAT).into());

        #[allow(unreachable_code)]
        {
            #[cfg(feature = "decimal")]
            if let Some(n) = rust_decimal::Decimal::from_f32(v) {
                return Ok(Dynamic::from_decimal(n));
            }

            Err(Error::custom(format!(
                "floating-point number is not supported: {v}"
            )))
        }
    }
    #[cfg(not(feature = "no_float"))]
    #[inline(always)]
    fn visit_f64<E: Error>(self, v: f64) -> Result<Self::Value, E> {
        #[cfg(not(feature = "no_float"))]
        return Ok((v as crate::FLOAT).into());

        #[allow(unreachable_code)]
        {
            #[cfg(feature = "decimal")]
            if let Some(n) = rust_decimal::Decimal::from_f64(v) {
                return Ok(Dynamic::from_decimal(n));
            }

            Err(Error::custom(format!(
                "floating-point number is not supported: {v}"
            )))
        }
    }

    #[cfg(feature = "no_float")]
    #[cfg(feature = "decimal")]
    #[inline]
    fn visit_f32<E: Error>(self, v: f32) -> Result<Self::Value, E> {
        use std::convert::TryFrom;

        rust_decimal::Decimal::try_from(v)
            .map(|v| v.into())
            .map_err(Error::custom)
    }
    #[cfg(feature = "no_float")]
    #[cfg(feature = "decimal")]
    #[inline]
    fn visit_f64<E: Error>(self, v: f64) -> Result<Self::Value, E> {
        use std::convert::TryFrom;

        rust_decimal::Decimal::try_from(v)
            .map(|v| v.into())
            .map_err(Error::custom)
    }

    #[inline(always)]
    fn visit_char<E: Error>(self, v: char) -> Result<Self::Value, E> {
        Ok(v.into())
    }
    #[inline(always)]
    fn visit_str<E: Error>(self, v: &str) -> Result<Self::Value, E> {
        Ok(v.into())
    }
    #[inline(always)]
    fn visit_string<E: Error>(self, v: String) -> Result<Self::Value, E> {
        Ok(v.into())
    }
    #[cfg(not(feature = "no_index"))]
    #[inline(always)]
    fn visit_bytes<E: Error>(self, v: &[u8]) -> Result<Self::Value, E> {
        Ok(Dynamic::from_blob(v.to_vec()))
    }

    #[inline(always)]
    fn visit_unit<E: Error>(self) -> Result<Self::Value, E> {
        Ok(Dynamic::UNIT)
    }

    #[inline(always)]
    fn visit_newtype_struct<D: Deserializer<'de>>(self, de: D) -> Result<Self::Value, D::Error> {
        Deserialize::deserialize(de)
    }

    #[cfg(not(feature = "no_index"))]
    fn visit_seq<A: serde::de::SeqAccess<'de>>(self, mut seq: A) -> Result<Self::Value, A::Error> {
        let mut arr = crate::Array::new();

        while let Some(v) = seq.next_element()? {
            arr.push(v);
        }

        Ok(arr.into())
    }

    #[cfg(not(feature = "no_object"))]
    fn visit_map<M: serde::de::MapAccess<'de>>(self, mut map: M) -> Result<Self::Value, M::Error> {
        let mut m = crate::Map::new();

        while let Some((k, v)) = map.next_entry()? {
            m.insert(k, v);
        }

        Ok(m.into())
    }
}

impl<'de> Deserialize<'de> for Dynamic {
    #[inline(always)]
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        deserializer.deserialize_any(DynamicVisitor)
    }
}

impl<'de> Deserialize<'de> for ImmutableString {
    #[inline(always)]
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        let s: String = Deserialize::deserialize(deserializer)?;
        Ok(s.into())
    }
}

impl<'de> Deserialize<'de> for Scope<'_> {
    #[inline(always)]
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        #[derive(Debug, Clone, Hash, Deserialize)]
        struct ScopeEntry {
            pub name: Identifier,
            pub value: Dynamic,
            #[serde(default)]
            pub is_constant: bool,
        }

        struct VecVisitor;

        impl<'de> Visitor<'de> for VecVisitor {
            type Value = Scope<'static>;

            fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
                formatter.write_str("a sequence")
            }

            #[inline]
            fn visit_seq<A>(self, mut access: A) -> Result<Self::Value, A::Error>
            where
                A: SeqAccess<'de>,
            {
                let mut scope = access
                    .size_hint()
                    .map_or_else(Scope::new, Scope::with_capacity);

                while let Some(ScopeEntry {
                    name,
                    value,
                    is_constant,
                }) = access.next_element()?
                {
                    if is_constant {
                        scope.push_constant_dynamic(name, value);
                    } else {
                        scope.push_dynamic(name, value);
                    }
                }

                Ok(scope)
            }
        }

        deserializer.deserialize_seq(VecVisitor)
    }
}
