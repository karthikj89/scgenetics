import glob
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import streamlit as st 
from streamlit_plotly_events import plotly_events
import pickle5 as pickle
from sklearn.preprocessing import MinMaxScaler
from urllib.request import urlretrieve

matrix_url = 'https://github.com/karthikj89/scgenetics/blob/master/www/scgwas_matrix.pkl?raw=true'
dst = 'scgwas_matrix.pkl'
urlretrieve(matrix_url, dst)

st.set_page_config(layout="centered")
st.title("Single Cell GWAS Explorer")
st.write("Size of circles indicate enrichment score, opacity of circles indicate -log(p-value).")

def get_trait_mapping(categories, trait_names, trait_name):
	trait_name = trait_name.split("_ldsc")[0]
	idx = None
	for i in range(len(trait_names)):
		if trait_name in trait_names[i]:
			idx = i
			break
	if idx is not None:
		return categories[idx]
	return "Other" 

def create_matrix():
	healthy_path = "./cell_type_programs"
	healthy_cell_proc_path = "./nmf"
	disease_path = "./disease_progression_programs"
	disease_cell_proc_path = "./nmf_disease_healthy"
	trait_cat_file = pd.read_csv("trait_category_mappings.csv")
	subcategories = trait_cat_file["Subcategory"].values
	trait_names = trait_cat_file["traits"].values

	tissues = []
	traits = []
	cell_types = []
	e_scores = []
	p_vals = []
	categories = []
	trait_categories = []
	enhancer_types = []

	for folder in glob.glob(healthy_path + "/*"):
		tissue_name = folder.split("/")[-1]
		print(folder)
		for cell_type_folder in glob.glob(folder + "/baselineLD_v2.1/*"):
			print(cell_type_folder)
			cell_type_name = cell_type_folder.split("/")[-1]
			for trait_file in glob.glob(cell_type_folder + "/*"):
				trait_name = trait_file.split("/")[-1].split("_postprocess.txt")[0]
				trait_cat = get_trait_mapping(subcategories, trait_names, trait_name)
				trait_matrix = pd.read_csv(trait_file, sep="\t")
				enhancer_names = trait_matrix.index.values
				trait_matrix = trait_matrix.values
				for i in range(len(trait_matrix)):
					enhancer_type = enhancer_names[i]
					e_score = float(trait_matrix[i, 3])
					p_value = float(trait_matrix[i, 5])
					tissues.append(tissue_name)
					traits.append(trait_name)
					cell_types.append(cell_type_name)
					categories.append("healthy")
					trait_categories.append(trait_cat)
					enhancer_types.append(enhancer_type)
					p_vals.append(p_value)
					e_scores.append(e_score)

	for folder in glob.glob(healthy_cell_proc_path + "/*"):
		tissue_name = folder.split("/")[-1]
		print(folder)
		for cell_type_folder in glob.glob(folder + "/baselineLD_v2.1/*"):
			print(cell_type_folder)
			cell_type_name = cell_type_folder.split("/")[-1]
			for trait_file in glob.glob(cell_type_folder + "/*"):
				trait_name = trait_file.split("/")[-1].split("_postprocess.txt")[0]
				trait_cat = get_trait_mapping(subcategories, trait_names, trait_name)
				trait_matrix = pd.read_csv(trait_file, sep="\t")
				enhancer_names = trait_matrix.index.values
				trait_matrix = trait_matrix.values
				for i in range(len(trait_matrix)):
					enhancer_type = enhancer_names[i]
					e_score = float(trait_matrix[i, 3])
					p_value = float(trait_matrix[i, 5])
					tissues.append(tissue_name)
					traits.append(trait_name)
					cell_types.append(cell_type_name)
					categories.append("healthy cellular process")
					trait_categories.append(trait_cat)
					enhancer_types.append(enhancer_type)
					p_vals.append(p_value)
					e_scores.append(e_score)

	for folder in glob.glob(disease_path + "/*"):
		tissue_name = folder.split("/")[-1]
		print(folder)
		for cell_type_folder in glob.glob(folder + "/baselineLD_v2.1/*"):
			print(cell_type_folder)
			cell_type_name = cell_type_folder.split("/")[-1]
			for trait_file in glob.glob(cell_type_folder + "/*"):
				trait_name = trait_file.split("/")[-1].split("_postprocess.txt")[0]
				trait_cat = get_trait_mapping(subcategories, trait_names, trait_name)
				try:
					trait_matrix = pd.read_csv(trait_file, sep="\t")
					enhancer_names = trait_matrix.index.values
					trait_matrix = trait_matrix.values
				except:
					continue
				for i in range(len(trait_matrix)):
					enhancer_type = enhancer_names[i]
					e_score = float(trait_matrix[i, 3])
					p_value = float(trait_matrix[i, 5])
					tissues.append(tissue_name)
					traits.append(trait_name)
					cell_types.append(cell_type_name)
					categories.append("disease")
					trait_categories.append(trait_cat)
					enhancer_types.append(enhancer_type)
					p_vals.append(p_value)
					e_scores.append(e_score)

	for folder in glob.glob(disease_cell_proc_path + "/*"):
		tissue_name = folder.split("/")[-1]
		print(folder)
		for cell_type_folder in glob.glob(folder + "/baselineLD_v2.1/*"):
			print(cell_type_folder)
			cell_type_name = cell_type_folder.split("/")[-1]
			for trait_file in glob.glob(cell_type_folder + "/*"):
				trait_name = trait_file.split("/")[-1].split("_postprocess.txt")[0]
				trait_cat = get_trait_mapping(subcategories, trait_names, trait_name)
				try:
					trait_matrix = pd.read_csv(trait_file, sep="\t")
					enhancer_names = trait_matrix.index.values
					trait_matrix = trait_matrix.values
				except:
					continue
				for i in range(len(trait_matrix)):
					enhancer_type = enhancer_names[i]
					e_score = float(trait_matrix[i, 3])
					p_value = float(trait_matrix[i, 5])
					tissues.append(tissue_name)
					traits.append(trait_name)
					cell_types.append(cell_type_name)
					categories.append("disease cellular process")
					trait_categories.append(trait_cat)
					enhancer_types.append(enhancer_type)
					p_vals.append(p_value)
					e_scores.append(e_score)


	data = {
		"tissue_category": categories,
		"trait_category": trait_categories,
		"tissue": tissues,
		"trait": traits,
		"cell_type": cell_types,
		"enhancer_type": enhancer_types,
		"e_score": e_scores,
		"p_value": p_vals
	}

	dataframe = pd.DataFrame.from_dict(data)
	dataframe.to_pickle("scgwas_matrix.pkl")

def generate_heatmap(scgwas_matrix, tissue_cat_selected, trait_cat_selected, tissue_selected, enhancer_type_selected):
	subset = scgwas_matrix[(scgwas_matrix.tissue_category==tissue_cat_selected)&(scgwas_matrix.tissue==tissue_selected)&(scgwas_matrix.enhancer_type==enhancer_type_selected)&(scgwas_matrix.trait_category==trait_cat_selected)]
	traits, celltypes = subset['trait'].values, subset['cell_type'].values
	traits = np.unique(traits)
	celltypes = np.unique(celltypes)
	
	trait_celltype_escores = np.zeros((len(traits), len(celltypes)))
	trait_celltype_pvals = np.zeros((len(traits), len(celltypes)))
	trait_celltype_pvals_transformed = np.zeros((len(traits), len(celltypes)))

	for i in range(len(traits)):
		trait_pairs = subset.loc[subset["trait"] == traits[i]]
		for j in range(len(celltypes)):
			if celltypes[j] in trait_pairs["cell_type"].values:
				trait_celltype_escores[i, j] = trait_pairs[trait_pairs["cell_type"] == celltypes[j]]["e_score"].values[0]
				trait_celltype_pvals[i, j] = trait_pairs[trait_pairs["cell_type"] == celltypes[j]]["p_value"].values[0]
				trait_celltype_pvals_transformed[i, j] = -np.log10(trait_pairs[trait_pairs["cell_type"] == celltypes[j]]["p_value"].values[0])
			else:
				trait_celltype_escores[i, j] = 0
				trait_celltype_pvals[i, j] = 0

	if trait_celltype_pvals_transformed.shape[0] > 0 and trait_celltype_pvals_transformed.shape[1] > 0:
		if trait_celltype_pvals_transformed.max() == trait_celltype_pvals_transformed.min():
			trait_celltype_pvals_transformed = (trait_celltype_pvals_transformed - trait_celltype_pvals_transformed.min()) / (1 - trait_celltype_pvals_transformed.min())
		else:
			trait_celltype_pvals_transformed = (trait_celltype_pvals_transformed - trait_celltype_pvals_transformed.min()) / (trait_celltype_pvals_transformed.max() - trait_celltype_pvals_transformed.min())

	x = []
	y = []
	escores = []
	pvals = []
	pvals_orig = []
	texts = []
	for i in range(len(traits)):
		for j in range(len(celltypes)):
			x.append(j)
			y.append(i)
			escore = trait_celltype_escores[i, j] * 2
			if escore < 0:
				escore = 0
			escores.append(escore)
			pvals.append(trait_celltype_pvals_transformed[i, j])
			pvals_orig.append(trait_celltype_pvals[i, j])
			texts.append("P-Value: " + str(trait_celltype_pvals[i, j]) + "<br>" + "E-Score: " + str(trait_celltype_escores[i, j]))

	fig = go.Figure(data=[go.Scatter(
	    x=x, y=y,
	    mode='markers',
	    text=texts,
	    marker=dict(
			color = ["#d62728"] * len(x) * len(y),
	    	opacity = pvals,
	    	size = escores
	    ),
	    hoverlabel=dict(bgcolor="white")
	)])

	fig.update_layout(
	    xaxis = dict(
	        tickmode = 'array',
			tickvals=[i for i in range(len(celltypes))], 
			ticktext=celltypes
	    )
	)

	fig.update_layout(
	    yaxis = dict(
	        tickmode = 'array',
			tickvals=[i for i in range(len(traits))], 
			ticktext=traits
	    )
	)

	fig.update_layout(
	    height=750,
	    width=750,
	)

	return fig

try:
	scgwas_matrix = pickle.load(open("scgwas_matrix.pkl","rb"))
except:
	urlretrieve(matrix_url, dst)
	scgwas_matrix = pickle.load(open("scgwas_matrix.pkl","rb"))
trait_categories = list(np.unique(scgwas_matrix["trait_category"].values))
tissue_categories = list(np.unique(scgwas_matrix["tissue_category"].values))
tissues = list(np.unique(scgwas_matrix["tissue"].values))
enhancer_types = list(np.unique(scgwas_matrix["enhancer_type"].values))

tissue_cat_selected = st.sidebar.selectbox("Select Tissue Category: ", tissue_categories, index=tissue_categories.index("healthy"))
trait_cat_selected = st.sidebar.selectbox("Select Trait Categories: ", trait_categories, index=trait_categories.index("Blood Biomarker - RBC"))
tissue_selected = st.sidebar.selectbox("Select Tissues: ", tissues, index=tissues.index("ICA_bonemarrow"))
enhancer_type_selected = st.sidebar.selectbox("Select Enhancer Type: ", enhancer_types, index=enhancer_types.index("ABC_Road_GI_BLD"))

if tissue_cat_selected and trait_cat_selected and tissue_selected and enhancer_type_selected:
	fig = generate_heatmap(scgwas_matrix, tissue_cat_selected, trait_cat_selected, tissue_selected, enhancer_type_selected)
	st.plotly_chart(fig)

st.markdown(
    """
    <style>
    [data-testid="stSidebar"][aria-expanded="true"] > div:first-child {
        width: 500px;
    }
    [data-testid="stSidebar"][aria-expanded="false"] > div:first-child {
        width: 500px;
        margin-left: -500px;
    }
    </style>
    """,
    unsafe_allow_html=True,
)

st.markdown(
        f"""
<style>
    .reportview-container .main .block-container{{
        margin-left: -275px;
    }}

</style>
""",
        unsafe_allow_html=True,
    )
