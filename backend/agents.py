import os
from crewai import Agent, LLM
from backend.tools import search_pubmed # Changed this import
from dotenv import load_dotenv

load_dotenv()

# 1. Define the Groq Brain
groq_llm = LLM(
    model="groq/llama-3.3-70b-versatile",
    api_key=os.getenv("GROQ_API_KEY"),
    temperature=0.2  
)

# 2. Define the Specialized Agents

classifier_agent = Agent(
    role='Clinical Symptom Classifier',
    goal='Identify and extract formal medical symptoms from natural language user input.',
    backstory="Expert clinical scribe converting patient descriptions into formal terms.",
    llm=groq_llm,
    allow_delegation=False,
    verbose=True
)

matcher_agent = Agent(
    role='Medical Research Researcher',
    goal='Search PubMed to find peer-reviewed literature.',
    backstory="Medical research scientist using PubMed for evidence-based info.",
    tools=[search_pubmed], # Clean, direct function call
    llm=groq_llm,
    allow_delegation=False,
    verbose=True
)

advisor_agent = Agent(
    role='Health Information Advisor',
    goal='Provide a clear summary with mandatory safety warnings.',
    backstory="Patient safety officer ensuring findings are non-diagnostic.",
    llm=groq_llm,
    allow_delegation=False,
    verbose=True
)