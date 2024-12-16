




function phy_par_degenrate()
    mN, mΛ, mΣ, mXi = 938.919, 1115.68, 1193.15, 1318.05
    mπ, mK = 138.039, 495.665
    meta = 547.51
    mch = [(mπ, mΣ), (mK, mN), (meta, mΛ), (mK, mXi) ]
    fϕ = 100.69
    pp = Dict(
        :mch => mch,
        :decons => fϕ,
        :threshold => [sum(mch[i]) for i in 1:4],
        :meson_mass => (mπ, mK)
    )
    return pp
end

function par_physical_channel()
    mp, mΣp, mΣ0, mn, mΣm, mΛ = [938.272, 1189.37, 1192.64, 939.565, 1197.45, 1115.68]
    mKm, mπm, mπ0, mKbar0, mπp = [493.68, 139.57, 134.977, 497.65, 139.570]
    mch = [(mKm, mp), (mπm, mΣp), (mπ0, mΣ0), (mKbar0, mn), (mπp, mΣm), (mπ0, mΛ)]
    fϕ = 100.69
    pp = Dict(
        :mch => mch,
        :decons => fϕ,
        :threshold => [sum(mch[i]) for i in 1:6],
        :meson_mass => (mπp, mπ0, mπm, mKbar0, mKm),
        :mB => [mp, mn, mΛ, mΣm, mΣ0, mΣp]
        )
    return pp
end

function par_physical_channel1()
    mN, mLam, mΣ, mXi = 938.919, 1115.68, 1193.15, 1318.05
    mπ, mK = 138.039, 495.665

    mp, mΣp, mΣ0, mn, mΣm, mΛ = [mN, mΣ, mΣ, mN, mΣ, mLam]
    mKm, mπm, mπ0, mKbar0, mπp = [mK, mπ, mπ, mK, mπ]
    mch = [(mKm, mp), (mπm, mΣp), (mπ0, mΣ0), (mKbar0, mn), (mπp, mΣm), (mπ0, mΛ)]
    fϕ = 100.69
    pp = Dict(
        :mch => mch,
        :decons => fϕ,
        :threshold => [sum(mch[i]) for i in 1:6],
        :meson_mass => (mπp, mπ0, mπm, mKbar0, mKm)
        )
    return pp
end
