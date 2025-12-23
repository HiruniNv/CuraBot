from crewai import Task

def create_tasks(user_input, classifier, matcher, advisor):
    t1 = Task(
        description=f"Extract symptoms from: '{user_input}'",
        expected_output="A list of medical symptoms.",
        agent=classifier
    )
    t2 = Task(
        description="Search PubMed for matching conditions.",
        expected_output="A summary of research matches.",
        agent=matcher,
        context=[t1]
    )
    t3 = Task(
        description="Synthesize findings into a safe report with a disclaimer.",
        expected_output="A final health insight report.",
        agent=advisor,
        context=[t2]
    )
    return [t1, t2, t3]