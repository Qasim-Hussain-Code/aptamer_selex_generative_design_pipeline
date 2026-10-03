"""Render every figure from retained tables; unavailable neural results remain absent."""
import argparse
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.metrics import precision_recall_curve
from common import ROOT, rows, table, settings


def save(fig, name):
    (ROOT / "figures").mkdir(exist_ok=True)
    for suffix in ["png", "pdf"]:
        fig.savefig(ROOT / f"figures/{name}.{suffix}", dpi=220, bbox_inches="tight")
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.parse_args()
    cfg = settings()
    plt.rcParams.update({"font.size":9,"axes.spines.top":False,"axes.spines.right":False})
    fig, ax = plt.subplots(figsize=(12,2.5))
    labels = ["Public rounds\nserial FASTQ", "Primer filters\nSQLite counts", "Assay partition\nedit families", "SELEX baselines\nMarkov null", "Canonical RNA\nstructure ablation", "Fixed test set\nfamily bootstrap"]
    for i,label in enumerate(labels):
        ax.text(i,0,label,ha="center",va="center",bbox=dict(boxstyle="round,pad=0.5",fc="#eaf0f4",ec="#54758a"))
        if i<5:
            ax.annotate("",(i+0.58,0),(i+0.42,0),arrowprops=dict(arrowstyle="->"))
    ax.text(3,-0.7,"RaptGen: hosted job prepared; no returned results. AptaDiff: licence blocker.",ha="center")
    ax.set(xlim=(-0.7,5.7),ylim=(-1,0.8)); ax.axis("off")
    save(fig,"01_pipeline")
    qc = list(rows(ROOT/"results/preprocessing_counts.tsv"))
    fig,axes=plt.subplots(1,2,figsize=(10,3.4))
    for ax,target in zip(axes,["tg2","integrin"]):
        data=sorted([r for r in qc if r["target"]==target],key=lambda r:int(r["round"]))
        x=[int(r["round"]) for r in data]
        for field,label in [("raw_reads","Raw"),("retained_reads","Retained"),("unique_sequences","Unique")]:
            ax.plot(x,[int(r[field]) for r in data],"o-",label=label)
        ax.set(title=target,xlabel="SELEX round",ylabel="Reads or unique sequences")
        ax.legend(frameon=False)
    fig.tight_layout();save(fig,"02_dataset_composition")
    metrics=list(rows(ROOT/"results/evaluation_metrics.tsv"))
    intervals=list(rows(ROOT/"results/bootstrap_intervals.tsv"))
    primary=[r for r in metrics if float(r["threshold"])==cfg["family_threshold"] and r["metric"]=="average_precision"]
    methods=["kmer_ridge","markov","kmer_logistic","frequency_observed","enrichment_observed"]
    fig,axes=plt.subplots(1,2,figsize=(11,3.8))
    for ax,target in zip(axes,["tg2","integrin"]):
        for offset,strategy,color in [(-0.2,"naive_sequence_random","#916a9b"),(0,"random","#ad8454"),(0.2,"family","#43768f")]:
            for i,method in enumerate(methods):
                point=next((r for r in primary if r["target"]==target and r["strategy"]==strategy and r["method"]==method),None)
                if point:
                    interval=next(r for r in intervals if r["target"]==target and r["strategy"]==strategy and r["method"]==method and float(r["threshold"])==cfg["family_threshold"])
                    value=float(point["value"])
                    legend={"naive_sequence_random":"Independent sequence split","random":"Matched exact-only exclusion","family":"Family exclusion"}
                    ax.plot(i+offset,value,"o",color=color,label=legend[strategy] if i==0 else None)
                    ax.vlines(i+offset,float(interval["lower_95"]),float(interval["upper_95"]),color=color)
        ax.set(xticks=range(len(methods)),xticklabels=["Enrichment\nridge","Markov","Assay\nlogistic","Observed\nfrequency","Observed\nenrichment"],ylim=(0,1.04),ylabel="Average precision",title=target)
        ax.legend(frameon=False)
    fig.tight_layout();save(fig,"03_random_vs_family")
    predictions=list(rows(ROOT/"results/held_out_predictions.tsv"))
    fig,axes=plt.subplots(1,2,figsize=(10,3.8))
    for ax,target in zip(axes,["tg2","integrin"]):
        for method in methods:
            data=[r for r in predictions if r["target"]==target and r["strategy"]=="family" and r["method"]==method and float(r["threshold"])==cfg["family_threshold"]]
            if data:
                precision,recall,_=precision_recall_curve([int(r["experimental_label"]) for r in data],[float(r["score"]) for r in data])
                ax.step(recall,precision,where="post",label=method)
        ax.set(xlabel="Recall",ylabel="Precision",title=target,ylim=(0,1.05)); ax.legend(fontsize=7,frameon=False)
    fig.tight_layout();save(fig,"04_experimental_discrimination")
    curves=list(rows(ROOT/"results/generation_budget_curves.tsv"))
    fig,axes=plt.subplots(1,2,figsize=(10,3.2))
    for ax,target in zip(axes,["tg2","integrin"]):
        for seed in cfg["generation_seeds"]:
            data=[r for r in curves if r["target"]==target and int(r["seed"])==seed]
            ax.plot([int(r["budget"]) for r in data],[int(r["positive_family_recovery"]) for r in data],"o-",label=str(seed))
        ax.set(title=target,xlabel="Unique Markov candidates",ylabel="Recovered positive test families",xscale="log")
        ax.legend(title="Generation seed",fontsize=7,frameon=False)
    fig.tight_layout();save(fig,"05_generation_budget")
    fig,axes=plt.subplots(1,2,figsize=(10,3.2))
    for ax,target in zip(axes,["tg2","integrin"]):
        # The plot uses a fixed first-500 sample from each seed. Full distributions are retained.
        from itertools import islice
        data=[r for seed in cfg["generation_seeds"] for r in islice(rows(ROOT/f"data/generated/{target}/markov_{seed}.novelty.tsv.gz"),500)]
        ax.scatter([float(r["maximum_identity"]) for r in data],[float(r["score"]) for r in data],s=3,alpha=0.2)
        ax.set(title=target,xlabel="Maximum edit identity to training",ylabel="Markov log probability per nucleotide")
    fig.tight_layout();save(fig,"06_novelty_ranking")
    fig,axes=plt.subplots(1,2,figsize=(9,3.5))
    for ax,target in zip(axes,["tg2","integrin"]):
        for i,method in enumerate(["kmer_ridge","kmer_structure_ridge","kmer_logistic","kmer_structure_logistic"]):
            interval=next((r for r in intervals if r["target"]==target and r["strategy"]=="family" and r["method"]==method and float(r["threshold"])==cfg["family_threshold"]),None)
            if interval:
                ax.plot(i,float(interval["estimate"]),"o",color="#43768f")
                ax.vlines(i,float(interval["lower_95"]),float(interval["upper_95"]),color="#43768f")
        ax.set(xticks=range(4),xticklabels=["SELEX\nsequence","SELEX\n+structure","Assay\nsequence","Assay\n+structure"],ylabel="Average precision",title=target,ylim=(0,1.04))
    fig.tight_layout();save(fig,"07_structure_ablation")
    resource=list(rows(ROOT/"logs/resource_usage.tsv"))
    by_stage={}
    for r in resource:
        previous=by_stage.get(r["stage"])
        record=dict(r)
        if previous:
            for field in ["peak_rss_bytes","peak_observed_disk_bytes"]:
                record[field]=max(int(previous[field]),int(r[field]))
        by_stage[r["stage"]]=record
    successful=list(by_stage.values())
    fig,axes=plt.subplots(1,3,figsize=(12,max(6,0.24*len(successful))))
    for ax,field,label,scale in zip(axes,["elapsed_seconds","peak_rss_bytes","peak_observed_disk_bytes"],["Elapsed seconds","Peak RSS MB","Project bytes GB"],[1,1e6,1e9]):
        ax.barh(range(len(successful)),[float(r[field])/scale for r in successful],color="#54758a")
        ax.set(yticks=range(len(successful)),yticklabels=[r["stage"] for r in successful] if field=="elapsed_seconds" else [],xlabel=label)
    fig.suptitle("Latest elapsed time per stage; maximum sampled RSS and disk across all attempts",fontsize=10)
    fig.tight_layout();save(fig,"08_resources")
    table(ROOT/"results/figure_manifest.tsv",[dict(figure=f"{i:02d}",script="scripts/16_figures.py",status="produced",source=source) for i,source in enumerate(["pipeline_configuration","preprocessing_counts.tsv","evaluation_metrics.tsv;bootstrap_intervals.tsv","held_out_predictions.tsv","generation_budget_curves.tsv","data/generated/*.novelty.tsv.gz","bootstrap_intervals.tsv;paired_differences.tsv","logs/resource_usage.tsv"],1)])
    print("Rendered eight figures; neural methods are absent because no artifacts have returned")


if __name__ == "__main__":
    main()
