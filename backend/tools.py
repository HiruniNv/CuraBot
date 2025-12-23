import os
from Bio import Entrez
from langchain.tools import tool
from dotenv import load_dotenv

load_dotenv()
Entrez.email = os.getenv("NCBI_EMAIL")

class PubMedTools:
    @tool("search_pubmed")
    def search_pubmed(query: str):
        """Searches PubMed for medical literature related to symptoms."""
        try:
            handle = Entrez.esearch(db="pubmed", term=query, retmax=3)
            record = Entrez.read(handle)
            handle.close()
            
            id_list = record["IdList"]
            if not id_list:
                return "No medical literature found for this query."

            handle = Entrez.efetch(db="pubmed", id=",".join(id_list), rettype="abstract", retmode="text")
            results = handle.read()
            handle.close()
            return results
        except Exception as e:
            return f"Error searching PubMed: {str(e)}"