    ## Efficiencies calculations
    # ==========================

    eta_cyclen = ((1+(1/lamb_ma1))*(h_3 - h_4) - (h_2 - h_1))/((1+(1/lamb_ma1))*h_3 - h_2); #Energetic efficiency of the cycle
    eta_mec = 1 - k_mec*((1+(1/lamb_ma1))*(h_3 - h_4) + (h_2 - h_1))/((1+(1/lamb_ma1))*(h_3 - h_4) - (h_2 - h_1)); #Mechanical efficiency
    eta_toten = eta_mec * eta_cyclen; #Total energetic efficiency

    eta_cyclex = ((1+(1/lamb_ma1))*(h_3 - h_4) - (h_2 - h_1))/((1+(1/lamb_ma1))*e_3 - e_2); #Exergetic efficiency of the cycle
    eta_rotex = ((1+(1/lamb_ma1))*(h_3 - h_4) - (h_2 - h_1))/((1+(1/lamb_ma1))*(e_3 - e_4) - (e_2 - e_1)); #Exergetic efficiency of the rotor assembly
    eta_combex = (1/f)*((1+(1/lamb_ma1))*e_3 - e_2)/((1+(1/lamb_ma1))*h_3 - h_2); #Exergetic efficiency of the combustion
    eta_totex = eta_mec * eta_cyclex * eta_combex; #Total exergetic efficiency

    ## Energetic losses
    # =================

    flow_air = P_e/(eta_mec*((1+(1/lamb_ma1))*(h_3 - h_4) - (h_2 - h_1)));  #Air mass flow rate [kg_air/s]

    loss_mec = k_mec*((1+(1/lamb_ma1))*(h_3 - h_4) + (h_2 - h_1))*flow_air; #Mechanical losses [W]
    loss_ech = ((1+(1/lamb_ma1))*h_4 - h_1)*flow_air;                       #Exhaust gas losses [W]

    ## Exergetic losses
    # =================

    loss_combex = flow_air*(e_2 + f*LHV*1e3/lamb_ma1 - (1+(1/lamb_ma1))*e_3); #Combustion losses [W]
    loss_turbine = flow_air*(1+(1/lamb_ma1))*((e_3 - e_4) - (h_3 - h_4)); #Shaft losses at the turbine [W]
    loss_compressor = flow_air*((h_2 - h_1) - (e_2 - e_1));               #Shaft losses at the compressor [W]
    loss_rotex = (loss_compressor + loss_turbine);                        #Total shaft losses [W]
    loss_echex = flow_air*((1+(1/lamb_ma1))*e_4 - e_1);                   #Exhaust gas losses [W]

    ## Mass flows
    # ===========

    dotm_a = flow_air; #Air mass flow rate [kg_air/s]
    dotm_f = flow_air/lamb_ma1; #Combustible mass flow rate [kg_comb/s]
    dotm_g = (1+(1/lamb_ma1))*flow_air;#Fluegas mass flow rate [kg_fluegas/s]

    gas_prop = [flue_conc_mass[2]*(1+(1/lamb_ma1))*flow_air,
                flue_conc_mass[3]*(1+(1/lamb_ma1))*flow_air,
                flue_conc_mass[1]*(1+(1/lamb_ma1))*flow_air,
                flue_conc_mass[0]*(1+(1/lamb_ma1))*flow_air]


    ## Generate graph to export:
    # ==========================

    # My 1st figure: Pie chart energy
    fig_pie_en = plt.figure(1)
    labels = 'Effective Power \n'+ '%.1f'%(P_e*1e-6)+' MW', 'Mechanical losses \n'+'%.1f'%(loss_mec*1e-6)+' MW', 'Exhaust losses MW \n'+'%.1f'%(loss_ech*1e-6)+' MW'
    sizes = [P_e*1e-6, loss_mec*1e-6, loss_ech*1e-6]
    colors = ['gold', 'yellowgreen', 'lightcoral']
    #explode = (0, 0, 0, 0)
    plt.pie(sizes, labels=labels, autopct='%1.1f%%', shadow=True, startangle=140)
    plt.axis('equal')
    plt.title("Primary energy flux " + "%.1f" %(LHV*dotm_f*1e-3) + " MW")

    # My 2nd figure: Pie chart exergy
    fig_pie_ex = plt.figure(2)
    labels = 'Effective Power \n'+ '%.1f'%(P_e*1e-6)+' MW', 'Mechanical losses \n'+'%.1f'%(loss_mec*1e-6)+' MW', 'Exhaust losses \n'+'%.1f'%(loss_echex*1e-6)+' MW', 'Turbine & compressor \n irreversibilities \n'+'%.1f'%(loss_rotex*1e-6)+' MW', 'Combustion \n irreversibilities \n'+'%.1f'%(loss_combex*1e-6)+' MW'
    sizes = [P_e*1e-6, loss_mec*1e-6, loss_echex*1e-6, loss_rotex*1e-6, loss_combex*1e-6]
    colors = ['gold', 'yellowgreen', 'lightcoral', 'lightskyblue', 'red']
    #explode = (0, 0, 0, 0)
    plt.pie(sizes, labels=labels, autopct='%1.1f%%', shadow=True, startangle=300)
    plt.axis('equal')
    plt.title("Primary exergy flux " + "%.1f" %(e_c*dotm_f*1e-3) + " MW")


    p = (p_1, p_2, p_3, p_4) # [Pa]
    T = (T_1, T_2, T_3, T_4) # [K]
    s = (s_1, s_2, s_3, s_4) # [J/kg/K]
    h = (h_1, h_2, h_3, h_4) # [J/kg]
    e = (e_1, e_2, e_3, e_4) # [J/kg]


    # My 3rd figure: T-s
    fig_Ts_diagram = plt.figure(3)
    plt.scatter(np.array(s)*1e-3, T, c="red")
    labels = ['1', '2', '3', '4'];
    for i, label in enumerate(labels): #get (0, label)
        plt.annotate(label,
                xy=(s[i]*1e-3, T[i]), #show point
                xytext=(5, 2), #show annotate
                textcoords='offset points',
                ha='right',
                va='bottom')
    plt.title('T-s diagram of the cycle')
    plt.xlabel("s $[kJ/kg/K]$")
    plt.ylabel("T $[K]$")


    # My 4th figure: h-s
    fig_hs_diagram = plt.figure(4)
    plt.scatter(np.array(s)*1e-3, np.array(h)*1e-3, c="red")
    labels = ['1', '2', '3', '4'];
    for i, label in enumerate(labels): #get (0, label)
        plt.annotate(label,
                xy=(s[i]*1e-3, h[i]*1e-3), #show point
                xytext=(5, 2), #show annotate
                textcoords='offset points',
                ha='right',
                va='bottom')
    plt.title('h-s diagram of the cycle')
    plt.xlabel("s $[kJ/kg/K]$")
    plt.ylabel("h $[kJ/kg]$")

    plt.show()
