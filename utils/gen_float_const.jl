p = BigFloat(pi)
e = exp(BigFloat(1))
table = [
    ["PI", p],
    ["E", e],
    ["LOG2_E", log2(e)],
    ["LN_10", log(big"10")],
    ["LN_2", log(big"2")],
    ["LOG10_2", log10(big"2")],
    ["LOG10_E", log10(e)],
    ["LOG2_10", log2(big"10")],
    ["FRAC_1_PI", 1 / p],
    ["FRAC_1_SQRT_2", 1 / sqrt(big"2")],
    ["FRAC_2_PI", 2 / p],
    ["FRAC_2_SQRT_PI", 2 / sqrt(p)],
    ["FRAC_PI_2", p / 2],
    ["FRAC_PI_3", p / 3],
    ["FRAC_PI_4", p / 4],
    ["FRAC_PI_6", p / 6],
    ["FRAC_PI_8", p / 8],
    ["SQRT_2", sqrt(big"2")],
    ["TAU", p * 2],
]

for (name, value) in table
    println("#[inline]")
    hi = Float64(value)
    lo = Float64(value - hi)
    println("fn $name() -> Self { Self { hi: $hi, lo: $lo } }")
    println("")
end
