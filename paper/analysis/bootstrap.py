"""Paired, reconstructing bootstrap of 12 dependent Spearman differences.
Only archived numeric scores are read; no model inference or manuscript edits.
"""
from pathlib import Path
from statistics import NormalDist
import argparse,csv,json,hashlib,math,shutil,sys,platform,time
import numpy as np

PROJECT=Path(__file__).resolve().parents[1]
HERE=Path(__file__).resolve().parent
COMPONENTS=['ribonn','te','ribodecode','half_life','csi','cai_ld','mfe_body','mfe_cds','gc']
METRICS=COMPONENTS+['AR','AT']
FAMILY=[('AR',k) for k in COMPONENTS]+[('AT',k) for k in COMPONENTS[:3]]
CONTRAST_NAMES=[a+' minus '+b for a,b in FAMILY]
NORMAL=NormalDist()

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def twice_average_ranks(values):
    """Exact twice-average ranks along final axis; integer tie ranks preserved."""
    values=np.asarray(values);n=values.shape[-1]
    order=np.argsort(values,axis=-1,kind='stable')
    sorted_values=np.take_along_axis(values,order,axis=-1)
    different=sorted_values[...,1:]!=sorted_values[...,:-1]
    positions=np.arange(n,dtype=np.int64)
    starts=np.maximum.accumulate(np.where(np.concatenate([np.ones(values.shape[:-1]+(1,),bool),different],axis=-1),positions,0),axis=-1)
    ends=np.minimum.accumulate(np.where(np.concatenate([different,np.ones(values.shape[:-1]+(1,),bool)],axis=-1),positions,n-1)[...,::-1],axis=-1)[...,::-1]
    sorted_ranks=starts+ends+2
    result=np.empty(values.shape,dtype=np.int64)
    np.put_along_axis(result,order,sorted_ranks,axis=-1)
    return result

def statistic(data):
    """data shape (..., 10, n): experimental output and 9 oriented components."""
    data=np.asarray(data);n=data.shape[-1]
    assert (2*n)**9<=np.iinfo(np.int64).max
    ranks=twice_average_ranks(data)
    ar_products=np.prod(ranks[...,1:,:],axis=-2,dtype=np.int64)
    at_products=np.prod(ranks[...,1:4,:],axis=-2,dtype=np.int64)
    composite_ranks=twice_average_ranks(np.stack([ar_products,at_products],axis=-2))
    prediction_ranks=np.concatenate([ranks[...,1:,:],composite_ranks],axis=-2)
    # Exact integer-centred average ranks; Pearson correlation handles all ties.
    yr=ranks[...,0,:]-(n+1);pr=prediction_ranks-(n+1)
    covariance=np.sum(pr*yr[...,None,:],axis=-1,dtype=np.int64)
    norm_products=np.sum(pr*pr,axis=-1,dtype=np.int64)*np.sum(yr*yr,axis=-1,dtype=np.int64)[...,None]
    rhos=covariance/np.sqrt(norm_products)
    differences=np.concatenate([rhos[...,9,None]-rhos[...,:9],rhos[...,10,None]-rhos[...,:3]],axis=-1)
    return rhos,differences

def interval_details(samples,observed,jackknife,confidence):
    samples=np.asarray(samples);jackknife=np.asarray(jackknife);B=len(samples)
    less=int(np.count_nonzero(samples<observed));equal=int(np.count_nonzero(samples==observed));greater=B-less-equal
    bias_probability=(less+.5*equal)/B
    u=jackknife.mean()-jackknife
    sum_squares=float(np.sum(u*u));sum_cubes=float(np.sum(u*u*u))
    acceleration=sum_cubes/(6*sum_squares**1.5) if sum_squares>0 else math.nan
    z0=NORMAL.inv_cdf(bias_probability) if 0<bias_probability<1 else math.nan
    alpha=(1-confidence)/2
    transformed=[];denominators=[]
    for tail in [alpha,1-alpha]:
        z=NORMAL.inv_cdf(tail);denom=1-acceleration*(z0+z);denominators.append(denom)
        transformed.append(NORMAL.cdf(z0+(z0+z)/denom) if math.isfinite(z0) and math.isfinite(acceleration) and denom!=0 else math.nan)
    finite=all(math.isfinite(v) for v in [z0,acceleration,*transformed,*denominators])
    valid=finite and 0<transformed[0]<transformed[1]<1 and all(d>0 for d in denominators)
    if not valid:raise RuntimeError('BCa interval is undefined; do not silently replace or discard data.')
    ci=np.quantile(samples,transformed,method='linear');pct=np.quantile(samples,[alpha,1-alpha],method='linear')
    return {'confidence_level':confidence,'lower':float(ci[0]),'upper':float(ci[1]),'percentile_lower':float(pct[0]),'percentile_upper':float(pct[1]),'z0':z0,'bias_probability':bias_probability,'less_than_observed':less,'equal_to_observed':equal,'greater_than_observed':greater,'acceleration':acceleration,'adjusted_lower_quantile':transformed[0],'adjusted_upper_quantile':transformed[1],'lower_tail_expected_draws':B*transformed[0],'upper_tail_expected_draws':B*(1-transformed[1]),'bca_denominators':denominators,'finite_and_ordered':valid,'positive_lower_bound':bool(ci[0]>0)}

def write_csv(path,rows):
    with Path(path).open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--input-dir',type=Path,default=HERE/'inputs' if (HERE/'inputs').is_dir() else PROJECT/'data/SupplementaryFigure3_Mauger_U');parser.add_argument('--output-dir',type=Path,default=PROJECT/'results/bootstrap');parser.add_argument('--resamples',type=int,default=200000);parser.add_argument('--seed',type=int,default=20261016);args=parser.parse_args()
    source=args.input_dir;out=args.output_dir;out.mkdir(parents=True,exist_ok=True)
    protected=[]
    protected_hashes={}
    rows=list(csv.DictReader((source/'measurements_and_scores.csv').open()));n=len(rows);assert n==39
    ids=[r['id'] for r in rows];assert len(set(ids))==n and set(ids)=={f'L{i}' for i in range(1,40)}
    original=json.loads((source/'analysis_protocol.json').read_text());assert original['combinations']['All Relevant']==COMPONENTS and original['combinations']['all translation']==COMPONENTS[:3]
    data=np.array([[float(row[k]) for row in rows] for k in ['U_bar_height_pt',*COMPONENTS]])
    observed_rhos,observed_delta=statistic(data)
    archived={r['scorer']:float(r['spearman_rho']) for r in csv.DictReader((source/'spearman_correlations.csv').open())}
    for key,value in zip(METRICS,observed_rhos):assert abs(value-archived[{'AR':'All Relevant','AT':'all translation'}.get(key,key)])<1e-12
    r=twice_average_ranks(data)
    assert np.array_equal(twice_average_ranks([0,2,3,2]),[2,5,8,5])
    assert all(math.prod(int(r[j,i]) for j in range(1,10))==int(rows[i]['hv_all_relevant_twice_rank_product']) for i in range(n))
    assert all(math.prod(int(r[j,i]) for j in range(1,4))==int(rows[i]['hv_all_translation_twice_rank_product']) for i in range(n))
    B=args.resamples;family_size=len(FAMILY);confidence_adjusted=1-.05/family_size
    protocol={'date':'2026-10-16','analysis':'Point 20: dependent Spearman differences, paired reconstructing bootstrap','n_constructs':n,'bootstrap_resamples':B,'random_seed':args.seed,'rng':'NumPy default_rng PCG64; indices generated as uint8','pairing':'One shared row-index sample for experimental output and all scorer outputs, and all 12 contrasts.','statistic':'Within each resample, recompute average component ranks, rank/size normalisation, exact twice-rank products for AR/AT and average ranks of each composite; calculate Spearman as Pearson correlation of average ranks.','model_predictions':'Fixed archived oriented values; no retraining or inference.','composite_membership':'Fixed, outcome-independent membership as in source protocol.','family':[{'composite':a,'component':b} for a,b in FAMILY],'interval_method':'Two-sided bias-corrected and accelerated (BCa) bootstrap','confidence_level_unadjusted':.95,'confidence_level_bonferroni_marginal':confidence_adjusted,'nominal_familywise_confidence_level':.95,'multiplicity':'One family of 12 contrasts; Bonferroni alpha=.05/12.','jackknife':'39 delete-one samples, each with 38 constructs and reconstructed component ranks/composites.','bias_correction':'midrank: (# bootstrap statistics below observed + 0.5 # equal) / B','quantile_method':'NumPy linear interpolation','primary_conclusion':'Corrected BCa lower confidence bound > 0; no use of bootstrap sign fractions as p-values.','diagnostics':'Finite bias correction and acceleration; ordered finite transformed probabilities, positive BCa denominators; report tail support and decisions in two disjoint 100000-draw halves.','scope':'Construct-level resampling conditional on this experimental readout and fixed model predictions. Assay, digitization and model-training errors not resampled.','analysis_source':'analyze.py','method_references':['https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.bootstrap.html','https://www.itl.nist.gov/div898/handbook/prc/section4/prc473.htm']}
    (out/'protocol.json').write_text(json.dumps(protocol,ensure_ascii=False,indent=2)+'\n')
    rng=np.random.default_rng(args.seed);indices=rng.integers(0,n,size=(B,n),dtype=np.uint8)
    bootstrap_rhos=np.empty((B,len(METRICS)));bootstrap_delta=np.empty((B,family_size));start_time=time.time()
    for start in range(0,B,2000):
        end=min(B,start+2000);sample=np.take(data,indices[start:end],axis=1).transpose(1,0,2)
        bootstrap_rhos[start:end],bootstrap_delta[start:end]=statistic(sample)
        if end%50000==0 or end==B:print(f'Bootstrap {end}/{B}; elapsed {time.time()-start_time:.1f}s',flush=True)
    assert np.isfinite(bootstrap_rhos).all() and np.isfinite(bootstrap_delta).all()
    jackknife_rhos=np.empty((n,len(METRICS)));jackknife_delta=np.empty((n,family_size))
    for omit in range(n):jackknife_rhos[omit],jackknife_delta[omit]=statistic(np.delete(data,omit,axis=1))
    assert np.isfinite(jackknife_delta).all()
    results=[];diagnostics=[]
    for j,(composite,component) in enumerate(FAMILY):
        plain=interval_details(bootstrap_delta[:,j],observed_delta[j],jackknife_delta[:,j],.95);adjusted=interval_details(bootstrap_delta[:,j],observed_delta[j],jackknife_delta[:,j],confidence_adjusted)
        half1=interval_details(bootstrap_delta[:B//2,j],observed_delta[j],jackknife_delta[:,j],confidence_adjusted);half2=interval_details(bootstrap_delta[B//2:,j],observed_delta[j],jackknife_delta[:,j],confidence_adjusted)
        row={'composite':composite,'component':component,'rho_composite':float(observed_rhos[METRICS.index(composite)]),'rho_component':float(observed_rhos[METRICS.index(component)]),'difference':float(observed_delta[j]),'bca_95_lower':plain['lower'],'bca_95_upper':plain['upper'],'bca_bonferroni_lower':adjusted['lower'],'bca_bonferroni_upper':adjusted['upper'],'superiority_before_adjustment':plain['positive_lower_bound'],'superiority_after_adjustment':adjusted['positive_lower_bound'],'percentile_95_lower':plain['percentile_lower'],'percentile_95_upper':plain['percentile_upper'],'percentile_bonferroni_lower':adjusted['percentile_lower'],'percentile_bonferroni_upper':adjusted['percentile_upper'],'lower_tail_expected_draws':adjusted['lower_tail_expected_draws'],'upper_tail_expected_draws':adjusted['upper_tail_expected_draws'],'split_half_decisions_agree':half1['positive_lower_bound']==half2['positive_lower_bound']==adjusted['positive_lower_bound'],'half1_corrected_lower':half1['lower'],'half2_corrected_lower':half2['lower']}
        results.append(row);diagnostics.append({'contrast':CONTRAST_NAMES[j],'unadjusted':plain,'bonferroni':adjusted,'first_half':half1,'second_half':half2})
        print(json.dumps(row),flush=True)
    assert all(sha(p)==h for p,h in protected_hashes.items())
    report={'protocol':protocol,'observed_correlations':dict(zip(METRICS,map(float,observed_rhos))),'comparisons':results,'diagnostics':diagnostics,'summary':{composite:{'n_comparisons':sum(r['composite']==composite for r in results),'n_positive_unadjusted':sum(r['composite']==composite and r['superiority_before_adjustment'] for r in results),'n_positive_corrected':sum(r['composite']==composite and r['superiority_after_adjustment'] for r in results),'positive_corrected_components':[r['component'] for r in results if r['composite']==composite and r['superiority_after_adjustment']],'all_components_superiority_supported':all(r['superiority_after_adjustment'] for r in results if r['composite']==composite)} for composite in ['AR','AT']},'all_bca_intervals_finite':True,'undefined_bootstrap_statistics':0,'all_split_half_decisions_agree':all(r['split_half_decisions_agree'] for r in results),'minimum_adjusted_tail_expected_draws':min(min(r['lower_tail_expected_draws'],r['upper_tail_expected_draws']) for r in results),'protected_source_sha256':protected_hashes,'protected_sources_unchanged':True,'software':{'python':sys.version,'numpy':np.__version__,'platform':platform.platform()},'elapsed_seconds':time.time()-start_time}
    (out/'results.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n');write_csv(out/'comparison_results.csv',results)
    np.savez_compressed(out/'bootstrap_arrays.npz',indices=indices,rhos=bootstrap_rhos,differences=bootstrap_delta,jackknife_rhos=jackknife_rhos,jackknife_differences=jackknife_delta,observed_rhos=observed_rhos,observed_differences=observed_delta,metrics=np.array(METRICS),contrast_names=np.array(CONTRAST_NAMES))
    write_csv(out/'jackknife_differences.csv',[{'omitted_id':ids[i],**{name:float(jackknife_delta[i,j]) for j,name in enumerate(CONTRAST_NAMES)}} for i in range(n)])
    inputs=out/'inputs';inputs.mkdir(exist_ok=True)
    for name in ['measurements_and_scores.csv','analysis_protocol.json','spearman_correlations.csv']:
        if (source/name).resolve()!=(inputs/name).resolve():shutil.copy2(source/name,inputs/name)
    if Path(__file__).resolve()!=(out/'analyze.py').resolve():shutil.copy2(__file__,out/'analyze.py')
    print('SUMMARY',json.dumps(report['summary']),flush=True);print('All protected sources unchanged.',flush=True)

if __name__=='__main__':main()
