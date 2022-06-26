    Ncomb=Decimal(composicoes[0])
    NO2=Decimal(composicoes[1])
    NH2O=Decimal(composicoes[2])
    NCO2=Decimal(composicoes[3])
    NCO=Decimal(composicoes[4])
    NO=Decimal(composicoes[5])
    NN=Decimal(composicoes[6])
    NNO2=Decimal(composicoes[7])
    NNO=Decimal(composicoes[8])
    NOH=Decimal(composicoes[9])
    NH2=Decimal(composicoes[10])
    NH=Decimal(composicoes[11])
    NN2 = Decimal(composicoes[12])
    #Net=composicoes[5]
    Ntot=Decimal(sum(composicoes))

    Ncomb=composicoes[0]
    NO2=composicoes[1]
    NH2O=composicoes[2]
    NCO2=composicoes[3]
    NCO=composicoes[4]
    NO=composicoes[5]
    NN=composicoes[6]
    NNO2=composicoes[7]
    NNO=composicoes[8]
    NOH=composicoes[9]
    NH2=composicoes[10]
    NH=composicoes[11]
    NN2=composicoes[12]
    #Net=composicoes[5]
    Ntot=sum(composicoes)


    deltaG=Decimal(Ncomb) * (gt_comb + (Decimal(Ncomb * P / Ntot).ln())) + \
           Decimal(NCO2) * (gt_co2 + (Decimal(NCO2 * P / Ntot).ln())) + \
           Decimal(NH2O) * (gt_h2o + (Decimal(NH2O * P / Ntot).ln())) +\
           Decimal(NO2) * (gt_o2 + (Decimal(NO2 * P / Ntot).ln()))+\
           Decimal(NO) * (gt_o + (Decimal(NO * P / Ntot).ln()))+\
           Decimal(NN) * (gt_n + (Decimal(NN * P / Ntot).ln()))+\
           Decimal(NNO2) * (gt_no2 + (Decimal(NNO2 * P / Ntot)))+\
           Decimal(NNO) * (gt_no + (Decimal(NNO * P / Ntot).ln()))+\
           Decimal(NOH) * (gt_oh + (Decimal(NOH * P / Ntot).ln()))+\
           Decimal(NH2) * (gt_h2 + (Decimal(NH2 * P / Ntot).ln()))+\
           Decimal(NH) * (gt_h + (Decimal(NH * P / Ntot).ln()))+\
           Decimal(NN2) * (gt_n2 + (Decimal(NN2 * P / Ntot).ln()))+\
           Decimal(NCO) * (gt_co + (Decimal(NCO * P / Ntot).ln()))

    gt_comb = Decimal(g0_comb+ float(sp.integrate(h_comb / (Ru * T ** 2), (T, 298.15, i))))
    #gt_gas = g0_gas - float(sp.integrate(h_gas / (Ru * T ** 2), (T, 298.15, i)))
    #gt_et = g0_et - float(sp.integrate(h_et / (Ru * T ** 2), (T, 298.15, i)))
    gt_co2 = Decimal(g0_co2 + float(sp.integrate(h_co2 / (Ru * T ** 2), (T, 298.15, i))))
    gt_o2 = Decimal(g0_o2 + float(sp.integrate(h_o2 / (Ru * T ** 2), (T, 298.15, i))))
    gt_co = Decimal(g0_co + float(sp.integrate(h_co / (Ru * T ** 2), (T, 298.15, i))))
    gt_o = Decimal(g0_o + float(sp.integrate(h_o / (Ru * T ** 2), (T, 298.15, i))))
    gt_n = Decimal(g0_n + float(sp.integrate(h_n / (Ru * T ** 2), (T, 298.15, i))))
    gt_no2 = Decimal(g0_no2 + float(sp.integrate(h_no2 / (Ru * T ** 2), (T, 298.15, i))))
    gt_no = Decimal(g0_no + float(sp.integrate(h_no / (Ru * T ** 2), (T, 298.15, i))))
    gt_oh = Decimal(g0_oh + float(sp.integrate(h_oh / (Ru * T ** 2), (T, 298.15, i))))
    gt_h2 = Decimal(g0_h2 + float(sp.integrate(h_h2 / (Ru * T ** 2), (T, 298.15, i))))
    gt_h = Decimal(g0_h + float(sp.integrate(h_h / (Ru * T ** 2), (T, 298.15, i))))
    gt_n2 = Decimal(g0_n2 + float(sp.integrate(h_n2 / (Ru * T ** 2), (T, 298.15, i))))

    deltaG=Ncomb * (gt_comb + float((Decimal(Ncomb) * Decimal(P) / Decimal(Ntot)).ln())) + \
           NCO2 * (gt_co2 + float((Decimal(NCO2) * Decimal(P) / Decimal(Ntot)).ln())) + \
           NH2O * (gt_h2o + float((Decimal(NH2O) * Decimal(P) / Decimal(Ntot)).ln())) +\
           NO2 * (gt_o2 + float((Decimal(NO2) * Decimal(P) / Decimal(Ntot)).ln()))+\
           NO * (gt_o + float((Decimal(NO) * Decimal(P) / Decimal(Ntot)).ln()))+\
           NN * (gt_n + float((Decimal(NN) * Decimal(P) / Decimal(Ntot)).ln()))+\
           NNO2 * (gt_no2 + float((Decimal(NNO2) * Decimal(P) / Decimal(Ntot)).ln()))+\
           NNO * (gt_no + float((Decimal(NNO) * Decimal(P) / Decimal(Ntot)).ln()))+\
           NOH * (gt_oh + float((Decimal(NOH) * Decimal(P) / Decimal(Ntot)).ln()))+\
           NH2 * (gt_h2 + float((Decimal(NH2) * Decimal(P) / Decimal(Ntot)).ln()))+\
           NH * (gt_h + float((Decimal(NH) * Decimal(P) / Decimal(Ntot)).ln()))+\
           NN2 * (gt_n2 + float((Decimal(NN2) * Decimal(P) / Decimal(Ntot)).ln()))+\
           NCO * (gt_co + float((Decimal(NCO) * Decimal(P) / Decimal(Ntot)).ln()))

    #print(deltaG)
    return abs(deltaG)

    deltaG=Ncomb * log(gt_comb + (Ncomb * P / Ntot)) + \
           NCO2 * (gt_co2 + log(NCO2 * P / Ntot)) + \
           NH2O * (gt_h2o + log(NH2O * P / Ntot)) +\
           NO2 * (gt_o2 + log(NO2 * P / Ntot))+\
           NO * (gt_o + log(NO * P / Ntot))+\
           NN * (gt_n + log(NN * P / Ntot))+\
           NNO2 * (gt_no2 + log(NNO2 * P / Ntot))+\
           NNO * (gt_no + log(NNO * P / Ntot))+\
           NOH * (gt_oh + log(NOH * P / Ntot))+\
           NH2 * (gt_h2 + log(NH2 * P / Ntot))+\
           NH * (gt_h + log(NH * P / Ntot))+\
           NN2 * (gt_n2 + log(NN2 * P / Ntot))+\
           NCO * (gt_co + log(NCO * P / Ntot))