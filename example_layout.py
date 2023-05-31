import plotly.graph_objects as go

LAYOUT = go.Layout(
    template="plotly_white",
    font=dict(
        color="black"
    ),
    legend=dict(
        font=dict(
            size=10
        )
    ),
    xaxis=dict(ticklen=0, linecolor="black"),
    yaxis=dict(ticklen=0, linecolor="black"),
    margin=dict(
            b=60,
            t=60,
    )
)