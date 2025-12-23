import os
from crewai import Agent, LLM
from .tools import PubMedTools
from dotenv import load_dotenv

# Load API keys from .env file
load_dotenv()

# 1. Define the Groq Brain (The LLM)
groq_llm = LLM(
    model="groq/llama-3.3-70b-versatile",
    api_key=os.getenv("GROQ_API_KEY"),
    temperature=0.2  
)

# 2. Define the Specialized Agents

# Agent 1: Extracts symptoms from messy user text
classifier_agent = Agent(
    role='Clinical Symptom Classifier',
    goal='Identify and extract formal medical symptoms from natural language user input.',
    backstory=(
        "You are an expert clinical scribe. Your job is to listen to patients and "
        "convert their descriptions (e.g., 'my head is pounding') into formal "
        "medical terms (e.g., 'severe throbbing headache')."
    ),
    llm=groq_llm,
    allow_delegation=False,
    verbose=True
)

# Agent 2: Searches PubMed for matches
matcher_agent = Agent(
    role='Medical Research Researcher',
    goal='Search PubMed to find peer-reviewed literature that matches the extracted symptoms.',
    backstory=(
        "You are a medical research scientist. You use the PubMed search tool to "
        "find evidence-based information. You look for potential conditions or "
        "medical explanations supported by recent studies."
    ),
    tools=[PubMedTools.search_pubmed],
    llm=groq_llm,
    allow_delegation=False,
    verbose=True
)

# Agent 3: Summarizes and adds safety disclaimers
advisor_agent = Agent(
    role='Health Information Advisor',
    goal='Provide a clear, patient-friendly summary of the findings with mandatory safety warnings.',
    backstory=(
        "You are a patient safety officer. Your priority is to ensure the user "
        "understands the research findings without thinking it is a final diagnosis. "
        "You always include a strong disclaimer advising the user to see a doctor."
    ),
    llm=groq_llm,
    allow_delegation=False,
    verbose=True
)