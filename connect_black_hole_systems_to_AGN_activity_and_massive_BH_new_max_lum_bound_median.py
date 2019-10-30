import sys
import h5py
sys.path.append('/home/aklantbhowmick/anaconda3/lib/python3.7/site-packages')
sys.path.append('/home/aklantbhowmick/anaconda3/lib/python3.7/site-packages/scalpy/')
sys.path.append('/home/aklantbhowmick/anaconda3/envs/nbodykit-env/lib/python3.6/site-packages/')
import mdot_to_Lbol
import arepo_package
import scipy.interpolate
import numpy
import os

def extract_ids(ids_to_be_extracted,complete_ids):
    #print(complete_ids)
    #print(ids_to_be_extracted)
    return numpy.where(complete_ids==ids_to_be_extracted)[0]
    
    
vec_extract_ids=numpy.vectorize(extract_ids,excluded=['complete_ids'])

#run='L25_n256'
#basePath='/ufrc/lblecha/aklantbhowmick/arepo_runs_aklant/'+run+'/output/'
#save_output_path='/ufrc/lblecha/aklantbhowmick/richness_output_AGN_activity/'

run='MBII'

simulation='MBII'

save_output_path='/ufrc/lblecha/aklantbhowmick/richness_output_AGN_activity/'
load_output_path_FOF='/ufrc/lblecha/aklantbhowmick/FOF_tag_outputs/'
load_output_path_snap='/ufrc/lblecha/aklantbhowmick/%s_BHs/'%simulation

if not os.path.exists(save_output_path):
        print("Making directory for storing richenss")
        os.makedirs(save_output_path)

radiative_efficiency=0.1
total_conv=mdot_to_Lbol.get_conversion_factor_arepo(radiative_efficiency)

redshift_space=[0.06,0.6,1.0,1.5,2.0,2.5,3.0,3.5,4]
log_LINKING_LENGTH_space=numpy.linspace(-2,2,40)
log_luminosity_cuts_space=[0,42.0,43.0,44.0,45.0]
log_bhmass_cuts_space=[6.0,7.0,8.0]

#log_lum_cut=.
#log_luminosity_cuts_space=[42.5,42.0]


#log_bhmass_cut=[6.]

for redshift in reversed(redshift_space):
    
    print('Starting redshift ',redshift)
    desired_redshift=redshift

    hf = h5py.File(load_output_path_snap+run+'bh_lum_host_halo_with_id_and_mass_z%.2f_and_sub_halo_id_added_velocity_all_BHs'%(redshift))

    hf_tag = h5py.File(load_output_path_snap+run+'bh_lum_host_halo_with_id_and_mass_z%.2f_and_sub_halo_id_added_velocity_all_BHs_tag_max_lum'%(redshift))

#    hf_half = h5py.File(load_output_path_snap+run+'bh_lum_host_halo_with_id_and_mass_z%.2f_and_sub_halo_id_added_velocity_all_BHs_division_eddington_ratio_max_lum_median'%(redshift))


    hf_half = h5py.File(load_output_path_snap+run+'bh_lum_host_halo_with_id_and_mass_z%.2f_and_sub_halo_id_added_velocity_all_BHs_division_max_lum_median'%(redshift))




    BH_Mass=hf.get('mass')[:]
    host_halo_mass=hf.get('host_halo_mass')[:]
    bolometric_luminosity=hf.get('L_bol')[:]
    host_subhalo_SM=hf.get('host_subhalo_SM')[:]

    ParticleIDs=hf.get('ids')[:]
    central_satellite_tag=hf_tag.get('central_satellite_tag')[:]

#    print(tag)

    half_by_halomass=hf_half.get('half_by_halomass')[:]
    half_by_bhmass=hf_half.get('half_by_bhmass')[:]
    half_by_stellarmass=hf_half.get('half_by_stellarmass')[:]
    


    
    for log_bhmass_cut in log_bhmass_cuts_space:
#        for log_LINKING_LENGTH in log_LINKING_LENGTH_space:
        scale_of_interest_space=numpy.array([0.015,0.1,1.0])
        for scale_of_interest in scale_of_interest_space:
        #for log_LINKING_LENGTH in log_LINKING_LENGTH_space:
                log_scale_of_interest=log_LINKING_LENGTH_space[log_LINKING_LENGTH_space<numpy.log10(scale_of_interest)][-1]
                log_LINKING_LENGTH=log_scale_of_interest
                final_FOF_tag,subsampled_quantities=numpy.load(load_output_path_FOF+run+'_log_bhmass_cut_%.1f_redshift_%.2f_log_LINKING_LENGTH_%.2f.npy'%(log_bhmass_cut,redshift,log_LINKING_LENGTH),allow_pickle=True)
                #print(final_FOF_tag)
                position_vector_cut,velocity_vector_gadget_units_cut,blackhole_id_cut,subhalo_id_cut,hosthalo_id_cut=subsampled_quantities
                complete_tag=numpy.array(list(zip(final_FOF_tag,hosthalo_id_cut)))

                unique_complete_tag=numpy.unique(complete_tag,axis=0)
                #unique_FOF_tag=numpy.unique(final_FOF_tag)
                RICHNESS=[]
                RICHNESS_BH_70=[]
                RICHNESS_BH_80=[]
                RICHNESS_BH_90=[]
                RICHNESS_BH_100=[]

                RICHNESS_active_42=[]
                RICHNESS_active_43=[]
                RICHNESS_active_44=[]
                RICHNESS_active_45=[]
                RICHNESS_active_46=[]
                RICHNESS_active_47=[]


                average_log_bolometric_luminosity=[]
                average_log_BH_Mass=[]

                primary_tag=[]
                primary_half_by_halomass=[]
                primary_half_by_bhmass=[]
                primary_half_by_stellarmass=[]
                BH_Mass_primary=[]
                host_halo_mass_primary=[]
                host_subhalo_SM_primary=[]

                bolometric_luminosity_primary=[]

                for tag,FOF_id in unique_complete_tag:
                    extract_FOFs=(final_FOF_tag==tag)&(FOF_id==hosthalo_id_cut)
                    ids_of_members=blackhole_id_cut[extract_FOFs]
                    bolometric_luminosity_members=bolometric_luminosity[vec_extract_ids(ids_to_be_extracted=ids_of_members,complete_ids=ParticleIDs)]
#                    extract_active_members=bolometric_luminosity_members>10**log_lum_cut                    


                    BH_Mass_members=BH_Mass[vec_extract_ids(ids_to_be_extracted=ids_of_members,complete_ids=ParticleIDs)]
                   
                    host_halo_mass_members=host_halo_mass[vec_extract_ids(ids_to_be_extracted=ids_of_members,complete_ids=ParticleIDs)]

                    host_subhalo_SM_members=host_subhalo_SM[vec_extract_ids(ids_to_be_extracted=ids_of_members,complete_ids=ParticleIDs)]

                    central_satellite_tag_members=central_satellite_tag[vec_extract_ids(ids_to_be_extracted=ids_of_members,complete_ids=ParticleIDs)]


                      
                    half_by_halomass_members=half_by_halomass[vec_extract_ids(ids_to_be_extracted=ids_of_members,complete_ids=ParticleIDs)]
                    
                    half_by_bhmass_members=half_by_bhmass[vec_extract_ids(ids_to_be_extracted=ids_of_members,complete_ids=ParticleIDs)]
 
                    half_by_stellarmass_members=half_by_stellarmass[vec_extract_ids(ids_to_be_extracted=ids_of_members,complete_ids=ParticleIDs)]

                    primary_bolometric_luminosity=numpy.max(bolometric_luminosity_members) 

                   
                    extract_primary=primary_bolometric_luminosity==bolometric_luminosity_members
                    bolometric_luminosity_primary.append(primary_bolometric_luminosity)
                    BH_Mass_primary.append((BH_Mass_members[extract_primary])[0])
                    host_halo_mass_primary.append((host_halo_mass_members[extract_primary])[0])
                    host_subhalo_SM_primary.append((host_subhalo_SM_members[extract_primary])[0])


                    primary_tag.append((central_satellite_tag_members[extract_primary])[0])
                    primary_half_by_halomass.append((half_by_halomass_members[extract_primary])[0])
                    primary_half_by_bhmass.append((half_by_bhmass_members[extract_primary])[0]) 
                    primary_half_by_stellarmass.append((half_by_stellarmass_members[extract_primary])[0])

                    RICHNESS.append(len(ids_of_members))

                    extract_active_members=BH_Mass_members>10**7
                    RICHNESS_BH_70.append(len(ids_of_members[extract_active_members]))

                    extract_active_members=BH_Mass_members>10**8
                    RICHNESS_BH_80.append(len(ids_of_members[extract_active_members]))

                    extract_active_members=BH_Mass_members>10**9
                    RICHNESS_BH_90.append(len(ids_of_members[extract_active_members]))

                    extract_active_members=BH_Mass_members>10**10
                    RICHNESS_BH_100.append(len(ids_of_members[extract_active_members]))

                    extract_active_members=bolometric_luminosity_members>10**42
                    RICHNESS_active_42.append(len(ids_of_members[extract_active_members]))

                    extract_active_members=bolometric_luminosity_members>10**43
                    RICHNESS_active_43.append(len(ids_of_members[extract_active_members]))

                    extract_active_members=bolometric_luminosity_members>10**44
                    RICHNESS_active_44.append(len(ids_of_members[extract_active_members]))

                    extract_active_members=bolometric_luminosity_members>10**45
                    RICHNESS_active_45.append(len(ids_of_members[extract_active_members]))

                    extract_active_members=bolometric_luminosity_members>10**46
                    RICHNESS_active_46.append(len(ids_of_members[extract_active_members]))

                    extract_active_members=bolometric_luminosity_members>10**47
                    RICHNESS_active_47.append(len(ids_of_members[extract_active_members]))


                                    
                RICHNESS=numpy.array(RICHNESS)
                primary_tag=numpy.array(primary_tag)
                primary_half_by_halomass=numpy.array(primary_half_by_halomass)   
               	
                primary_half_by_bhmass=numpy.array(primary_half_by_bhmass)   
                BH_Mass_primary=numpy.array(BH_Mass_primary)               	
                bolometric_luminosity_primary=numpy.array(bolometric_luminosity_primary)
                host_halo_mass_primary=numpy.array(host_halo_mass_primary)

                RICHNESS_BH_70=numpy.array(RICHNESS_BH_70)
                RICHNESS_BH_80=numpy.array(RICHNESS_BH_80)
                RICHNESS_BH_90=numpy.array(RICHNESS_BH_90)
                RICHNESS_BH_100=numpy.array(RICHNESS_BH_100)


                RICHNESS_active_42=numpy.array(RICHNESS_active_42)
                RICHNESS_active_43=numpy.array(RICHNESS_active_43)
                RICHNESS_active_44=numpy.array(RICHNESS_active_44)
                RICHNESS_active_45=numpy.array(RICHNESS_active_45)
                RICHNESS_active_46=numpy.array(RICHNESS_active_46)
                RICHNESS_active_47=numpy.array(RICHNESS_active_47)


                #average_log_bolometric_luminosity=numpy.array(average_log_bolometric_luminosity)
                #average_log_BH_Mass=numpy.array(average_log_BH_Mass)
                
#                numpy.save(save_output_path+run+'log_bhmass_cut_%.1f_redshift_%.2f_log_LINKING_LENGTH_%.2f_all_info_more_split_by_Ledd_bound_median.npy'%(log_bhmass_cut,redshift,log_LINKING_LENGTH),[RICHNESS,RICHNESS_active_42,RICHNESS_active_43,RICHNESS_active_44,RICHNESS_active_45,RICHNESS_active_46,RICHNESS_active_47,RICHNESS_BH_70,RICHNESS_BH_80,RICHNESS_BH_90,RICHNESS_BH_100,primary_tag,primary_half_by_halomass,primary_half_by_bhmass,primary_half_by_stellarmass,BH_Mass_primary,bolometric_luminosity_primary,host_halo_mass_primary,host_subhalo_SM_primary])


                numpy.save(save_output_path+run+'log_bhmass_cut_%.1f_redshift_%.2f_log_LINKING_LENGTH_%.2f_all_info_more_bound_median.npy'%(log_bhmass_cut,redshift,log_LINKING_LENGTH),[RICHNESS,RICHNESS_active_42,RICHNESS_active_43,RICHNESS_active_44,RICHNESS_active_45,RICHNESS_active_46,RICHNESS_active_47,RICHNESS_BH_70,RICHNESS_BH_80,RICHNESS_BH_90,RICHNESS_BH_100,primary_tag,primary_half_by_halomass,primary_half_by_bhmass,primary_half_by_stellarmass,BH_Mass_primary,bolometric_luminosity_primary,host_halo_mass_primary,host_subhalo_SM_primary])
 
        #print(RICHNESS)
        #print(RICHNESS_active)
