static  void contact_point_propagate(
        Front           *front,
        POINTER         wave,
        POINT           *oldp,
        POINT           *newp,
        HYPER_SURF_ELEMENT  *oldhse,
        HYPER_SURF      *oldhs,
        double          dt,
        double          *V)
{
    int i, dim = front->rect_grid->dim;
    double **m_mom = eqn_params->mom;
    double *m_dens = eqn_params->dens;
    double *m_engy = eqn_params->engy;
    double *m_pres = eqn_params->pres;
    double *p0;
    STATE *oldst,*newst;
    STATE *newstl, *newstr;
    POINTER sl,sr;
    EOS_PARAMS *eos = eqn_params->eos;
    double default_var;

    // JAM - Don't need this
    //FT_GetStatesAtPoint(oldp,oldhse,oldhs,&sl,&sr);

        newstl = (STATE*)left_state(newp);
        newstl->eos = &eos[negative_component(oldhs)];

        newstr = (STATE*)right_state(newp);
        newstr->eos = &eos[positive_component(oldhs)];

    STATE *stl, *str;
    STATE *ssl, *ssr;
    FT_ScalarMemoryAlloc((POINTER*)&stl,sizeof(STATE));
    FT_ScalarMemoryAlloc((POINTER*)&str,sizeof(STATE));
    FT_ScalarMemoryAlloc((POINTER*)&ssl,sizeof(STATE));
    FT_ScalarMemoryAlloc((POINTER*)&ssr,sizeof(STATE));
        *stl = *(STATE*)sl;
        *str = *(STATE*)sr;
        *ssl = *(STATE*)sl;
        *ssr = *(STATE*)sr;
        stl->eos = newstl->eos;
        str->eos = newstr->eos;

    //TODO: Replace this with the stuff from FLASH grid
    if (!g_intrp_state(stl,str,oldp,oldhs,front))
    {
        //if can't find block for the point, keep the old point
        *newp = *oldp;
        return;
    }

        double lvnor, rvnor;
        double lvtang[MAXD], rvtang[MAXD];
        double nor[MAXD];
        double velp[MAXD],s;
        const double epsilon = 1.e-10;
        const double delta = 1.e-10;

        // JAM - Do I need this?
        normal(oldp,oldhse,oldhs,nor,front);

        lvnor = 0;
        rvnor = 0;
        for (i=0; i<dim; i++)
        {
                lvnor += stl->vel[i]*nor[i];
                rvnor += str->vel[i]*nor[i];
        }
        for (i=0; i<dim; i++)
        {
                lvtang[i] = stl->vel[i]-lvnor*nor[i];
                rvtang[i] = str->vel[i]-rvnor*nor[i];
        }

    double pml, pmr, uml, umr, ml, mr;
    RIEMANN_SOLVER_WAVE_TYPE l_wave,r_wave;
    STATE *stl2, *str2;
    STATE *ansl, *ansr;

    FT_ScalarMemoryAlloc((POINTER*)&stl2,sizeof(STATE));
    FT_ScalarMemoryAlloc((POINTER*)&str2,sizeof(STATE));
    *stl2 = *stl;
    *str2 = *str;
        stl2->eos = newstl->eos;
        str2->eos = newstr->eos;

    stl2->vel[0] = lvnor;
    str2->vel[0] = rvnor;
    for (i = 1; i < dim; i++)
        stl2->vel[i] = str2->vel[i] = 0.0;

    FT_ScalarMemoryAlloc((POINTER*)&ansl,sizeof(STATE));
    FT_ScalarMemoryAlloc((POINTER*)&ansr,sizeof(STATE));

    if (!find_mid_state(stl2,str2,0.0/*pjump*/,&pml,&pmr,&uml,&umr,&ml,&mr,
            &l_wave,&r_wave))
    {
        printf("In contact_point_propagate(), find_mid_state() failed at (%lf, %lf, %lf).\n",
            Coords(oldp)[0], Coords(oldp)[1], Coords(oldp)[2]);
        printf("nor: (%lf, %lf, %lf).\n", nor[0], nor[1], nor[2]);
        printf("ssl: dens = %e, pres = %e, vel = %e.\n",
            ssl->dens, ssl->pres, ssl->vel[0]);
        printf("ssr: dens = %e, pres = %e, vel = %e.\n",
            ssr->dens, ssr->pres, ssr->vel[0]);
        printf("stl: dens = %e, pres = %e, vel = %e.\n",
            stl2->dens, stl2->pres, stl2->vel[0]);
        printf("str: dens = %e, pres = %e, vel = %e.\n",
            str2->dens, str2->pres, str2->vel[0]);
        clean_up(ERROR);
    }
    midstate(stl2,ansl,ml,uml,pml,l_wave,1);
    midstate(str2,ansr,mr,umr,pmr,r_wave,-1);   //1 for left, -1 for right

        for (i=0; i<dim; i++)
            velp[i] = ansl->vel[0]*nor[i] + lvtang[i];

        newstl->dens = ansl->dens;
        newstl->pres = ansl->pres;
        newstr->dens = ansr->dens;
        newstr->pres = ansr->pres;
        if(eqn_params->multi_comp_non_reactive == YES)
        {
            int ii;
            for(ii = 0; ii < eqn_params->n_comps; ii++)
            {
                newstl->pdens[ii] = ansl->pdens[ii];
                newstr->pdens[ii] = ansr->pdens[ii];
            }
        }
        for (i=0; i<dim; i++)
        {
            newstl->vel[i] = velp[i];
            newstl->momn[i] = newstl->vel[i]*newstl->dens;
            newstr->vel[i] = velp[i];
            newstr->momn[i] = newstr->vel[i]*newstr->dens;
        }
        newstl->engy = EosEnergy(newstl);
        newstr->engy = EosEnergy(newstr);

    FT_FreeThese(2,stl2,str2);
    FT_FreeThese(2,ansl,ansr);

        for (i = 0; i < dim; ++i)
        {
            Coords(newp)[i] = Coords(oldp)[i] + dt*velp[i];
            set_max_front_speed(i,fabs(velp[i]),NULL,Coords(newp),front);
        }
        s = mag_vector(velp,dim);

    double cl, cr, c;
    cl = EosSoundSpeed(newstl);
    cr = EosSoundSpeed(newstr);
    c = std::max(cl, cr);
    set_max_front_speed(dim,s+c,NULL,Coords(newp),front);
    FT_FreeThese(4,stl,str,ssl,ssr);
}   /* end contact_point_propagate */

