import os
from Bio import Entrez
from crewai.tools import tool 
from dotenv import load_dotenv

load_dotenv()
Entrez.email = os.getenv("NCBI_EMAIL", "your.email@example.com")

@tool("search_pubmed")
def search_pubmed(query: str):
    """
    Searches PubMed for medical literature related to symptoms. 
    The query should be clinical keywords (e.g., 'acute tachycardia human trials').
    """
    try:
        # Refine query to avoid irrelevant results
        refined_query = f"{query} AND (clinical trial[Filter] OR review[Filter])"
        handle = Entrez.esearch(db="pubmed", term=refined_query, retmax=3)
        record = Entrez.read(handle)
        handle.close()
        
        id_list = record["IdList"]
        if not id_list:
            return "No high-quality medical literature found for this specific query."

        handle = Entrez.efetch(db="pubmed", id=",".join(id_list), rettype="abstract", retmode="text")
        results = handle.read()
        handle.close()
        return results
    except Exception as e:
        return f"Error searching PubMed: {str(e)}"
   