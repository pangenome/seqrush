op ==(int, int) -> bool;
op !=(int, int) -> bool;
op >(int, int) -> bool;
op >=(int, int) -> bool;
op <(int, int) -> bool;
op <=(int, int) -> bool;
op &(int, int) -> int;
op |(int, int) -> int;
op ^(int, int) -> int;
op ..(int, int) -> Range<int>;
op ..=(int, int) -> RangeInclusive<int>;

op ==(bool, bool) -> bool;
op !=(bool, bool) -> bool;
op >(bool, bool) -> bool;
op >=(bool, bool) -> bool;
op <(bool, bool) -> bool;
op <=(bool, bool) -> bool;
op &(bool, bool) -> bool;
op |(bool, bool) -> bool;
op ^(bool, bool) -> bool;

op ==((), ()) -> bool;
op !=((), ()) -> bool;
op >((), ()) -> bool;
op >=((), ()) -> bool;
op <((), ()) -> bool;
op <=((), ()) -> bool;

op +(int, int) -> int;
op -(int, int) -> int;
op *(int, int) -> int;
op /(int, int) -> int;
op %(int, int) -> int;
op **(int, int) -> int;
op >>(int, int) -> int;
op <<(int, int) -> int;

op +(float, float) -> float;
op -(float, float) -> float;
op *(float, float) -> float;
op /(float, float) -> float;
op %(float, float) -> float;
op **(float, float) -> float;
op ==(float, float) -> bool;
op !=(float, float) -> bool;
op >(float, float) -> bool;
op >=(float, float) -> bool;
op <(float, float) -> bool;
op <=(float, float) -> bool;

op +(float, int) -> float;
op -(float, int) -> float;
op *(float, int) -> float;
op /(float, int) -> float;
op %(float, int) -> float;
op **(float, int) -> float;
op ==(float, int) -> bool;
op !=(float, int) -> bool;
op >(float, int) -> bool;
op >=(float, int) -> bool;
op <(float, int) -> bool;
op <=(float, int) -> bool;

op +(int, float) -> float;
op -(int, float) -> float;
op *(int, float) -> float;
op /(int, float) -> float;
op %(int, float) -> float;
op **(int, float) -> float;
op ==(int, float) -> bool;
op !=(int, float) -> bool;
op >(int, float) -> bool;
op >=(int, float) -> bool;
op <(int, float) -> bool;
op <=(int, float) -> bool;

op +(Decimal, Decimal) -> Decimal;
op -(Decimal, Decimal) -> Decimal;
op *(Decimal, Decimal) -> Decimal;
op /(Decimal, Decimal) -> Decimal;
op %(Decimal, Decimal) -> Decimal;
op **(Decimal, Decimal) -> Decimal;
op ==(Decimal, Decimal) -> bool;
op !=(Decimal, Decimal) -> bool;
op >(Decimal, Decimal) -> bool;
op >=(Decimal, Decimal) -> bool;
op <(Decimal, Decimal) -> bool;
op <=(Decimal, Decimal) -> bool;

op +(Decimal, int) -> Decimal;
op -(Decimal, int) -> Decimal;
op *(Decimal, int) -> Decimal;
op /(Decimal, int) -> Decimal;
op %(Decimal, int) -> Decimal;
op **(Decimal, int) -> Decimal;
op ==(Decimal, int) -> bool;
op !=(Decimal, int) -> bool;
op >(Decimal, int) -> bool;
op >=(Decimal, int) -> bool;
op <(Decimal, int) -> bool;
op <=(Decimal, int) -> bool;

op +(int, Decimal) -> Decimal;
op -(int, Decimal) -> Decimal;
op *(int, Decimal) -> Decimal;
op /(int, Decimal) -> Decimal;
op %(int, Decimal) -> Decimal;
op **(int, Decimal) -> Decimal;
op ==(int, Decimal) -> bool;
op !=(int, Decimal) -> bool;
op >(int, Decimal) -> bool;
op >=(int, Decimal) -> bool;
op <(int, Decimal) -> bool;
op <=(int, Decimal) -> bool;

op +(String, String) -> String;
op -(String, String) -> String;
op ==(String, String) -> bool;
op !=(String, String) -> bool;
op >(String, String) -> bool;
op >=(String, String) -> bool;
op <(String, String) -> bool;
op <=(String, String) -> bool;

op +(char, char) -> String;
op ==(char, char) -> bool;
op !=(char, char) -> bool;
op >(char, char) -> bool;
op >=(char, char) -> bool;
op <(char, char) -> bool;
op <=(char, char) -> bool;

op +(char, String) -> String;
op ==(char, String) -> bool;
op !=(char, String) -> bool;
op >(char, String) -> bool;
op >=(char, String) -> bool;
op <(char, String) -> bool;
op <=(char, String) -> bool;

op +(String, char) -> String;
op -(String, char) -> String;
op ==(String, char) -> bool;
op !=(String, char) -> bool;
op >(String, char) -> bool;
op >=(String, char) -> bool;
op <(String, char) -> bool;
op <=(String, char) -> bool;

op +((), String) -> String;
op ==((), String) -> bool;
op !=((), String) -> bool;
op >((), String) -> bool;
op >=((), String) -> bool;
op <((), String) -> bool;
op <=((), String) -> bool;

op +(String, ()) -> String;
op ==(String, ()) -> bool;
op !=(String, ()) -> bool;
op >(String, ()) -> bool;
op >=(String, ()) -> bool;
op <(String, ()) -> bool;
op <=(String, ()) -> bool;

op +(Blob, Blob) -> Blob;
op +(Blob, char) -> Blob;
op ==(Blob, Blob) -> bool;
op !=(Blob, Blob) -> bool;


op ==(Range<int>, RangeInclusive<int>) -> bool;
op !=(Range<int>, RangeInclusive<int>) -> bool;

op ==(RangeInclusive<int>, Range<int>) -> bool;
op !=(RangeInclusive<int>, Range<int>) -> bool;

op ==(Range<int>, Range<int>) -> bool;
op !=(Range<int>, Range<int>) -> bool;

op ==(RangeInclusive<int>, RangeInclusive<int>) -> bool;
op !=(RangeInclusive<int>, RangeInclusive<int>) -> bool;

op ==(?, ?) -> bool;
op !=(?, ?) -> bool;
op >(?, ?) -> bool;
op >=(?, ?) -> bool;
op <(?, ?) -> bool;
op <=(?, ?) -> bool;


op &=(bool, bool);
op |=(bool, bool);

op +=(int, int);
op -=(int, int);
op *=(int, int);
op /=(int, int);
op %=(int, int);
op **=(int, int);
op >>=(int, int);
op <<=(int, int);
op &=(int, int);
op |=(int, int);
op ^=(int, int);

op +=(float, float);
op -=(float, float);
op *=(float, float);
op /=(float, float);
op %=(float, float);
op **=(float, float);

op +=(float, int);
op -=(float, int);
op *=(float, int);
op /=(float, int);
op %=(float, int);
op **=(float, int);

op +=(Decimal, Decimal);
op -=(Decimal, Decimal);
op *=(Decimal, Decimal);
op /=(Decimal, Decimal);
op %=(Decimal, Decimal);
op **=(Decimal, Decimal);

op +=(Decimal, int);
op -=(Decimal, int);
op *=(Decimal, int);
op /=(Decimal, int);
op %=(Decimal, int);
op **=(Decimal, int);

op +=(String, String);
op -=(String, String);
op +=(String, char);
op -=(String, char);
op +=(char, String);
op +=(char, char);

op +=(Array, Array);
op +=(Array, ?);

op +=(Blob, Blob);
op +=(Blob, int);
op +=(Blob, char);
op +=(Blob, String);

op in(?, Array) -> bool;
op in(String, String) -> bool;
op in(char, String) -> bool;
op in(int, Range<int>) -> bool;
op in(int, RangeInclusive<int>) -> bool;
op in(String, Map) -> bool;
op in(int, Blob) -> bool;
