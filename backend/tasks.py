from crewai import Task

def create_tasks(user_input, classifier, matcher, advisor):
    # Task 1: Refine symptoms into search queries
    t1 = Task(
        description=(
            f"Analyze the user's input: '{user_input}'. "
            "Identify the most likely clinical conditions and search terms."
        ),
        expected_output="A list of professional medical search keywords.",
        agent=classifier
    )
    
    # Task 2: Search for relevant clinical evidence
    t2 = Task(
        description=(
            "Using the refined keywords, search PubMed for recent clinical studies or trials. "
            "Focus on papers that explain the most likely causes of the symptoms."
        ),
        expected_output="A summary of the most relevant research findings.",
        agent=matcher,
        context=[t1]
    )
    
    # Task 3: Format the final output for the CuraBot Dashboard
    t3 = Task(
        description=(
            "Synthesize findings into a professional report using this exact structure:\n"
            "### 🩺 CURABOT MEDICAL RESEARCH SUMMARY\n"
            "Reported Symptoms: [Summary]\n\n"
            "--- \n"
            "### 🔍 KEY RESEARCH FINDINGS\n"
            "[Relevant study summaries]\n\n"
            "--- \n"
            "### 💡 CLINICAL CONSIDERATIONS\n"
            "[Contextualized medical reasoning]\n\n"
            "### 📋 NEXT STEPS\n"
            "[Professional advice and disclaimer]"
        ),
        expected_output="A complete Markdown-formatted medical research report.",
        agent=advisor,
        context=[t2]
    )
    
    return [t1, t2, t3]
