from crewai import Task

def create_tasks(user_input, classifier, matcher, advisor):
    t1 = Task(
        description=(
            f"Analyze the user's input: '{user_input}'. "
            "Identify and extract specific medical symptoms and duration. "
            "Convert them into professional medical keywords for a PubMed search."
        ),
        expected_output="A structured list of medical symptoms and search keywords.",
        agent=classifier
    )
    
    t2 = Task(
        description=(
            "Using the keywords from the classifier, search PubMed for recent clinical studies, "
            "systematic reviews, or trials. Focus strictly on human health data. "
            "Avoid results that are purely geographical or unrelated to the pathology of the symptoms."
        ),
        expected_output="A summarized list of the most relevant research papers with titles and key findings.",
        agent=matcher,
        context=[t1]
    )
    
    t3 = Task(
        description=(
            "Synthesize the research findings into a professional, empathetic report. "
            "Use the following structure strictly:\n"
            "### 🩺 CURABOT MEDICAL RESEARCH SUMMARY\n"
            "**Reported Symptoms:** [Summary of symptoms]\n\n"
            "--- \n"
            "### 🔍 KEY RESEARCH FINDINGS\n"
            "[List specific study insights here]\n\n"
            "--- \n"
            "### 💡 CLINICAL CONSIDERATIONS\n"
            "[General medical context based on the research]\n\n"
            "### 📋 NEXT STEPS\n"
            "- [Non-prescriptive advice like 'Consult a specialist']\n\n"
            "**Disclaimer:** CuraBot is an AI research tool. Not a diagnosis. Consult a doctor."
        ),
        expected_output="A complete, professional Markdown-formatted medical research report.",
        agent=advisor,
        context=[t2]
    )
    
    return [t1, t2, t3]
