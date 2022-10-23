using Dash
using DashBootstrapComponents

# https://dash-bootstrap-components.opensource.faculty.ai/docs/components/alert/
app = dash(external_stylesheets=[dbc_themes.BOOTSTRAP])

app.layout = dbc_container(
    dbc_alert("Hello boostrap", color="success"),
    className="p-5"
)
# app.layout = html_div() do
#     html_h1("Hello Dash"),
#     html_div("Dash: A web application framework for Julia"),
#     dcc_graph(
#         id = "example-graph-1",
#         figure = (
#             data = [
#                 (x = ["giraffes", "orangutans", "monkeys"], y = [20, 14, 23], type = "bar", name = "SF"),
#                 (x = ["giraffes", "orangutans", "monkeys"], y = [12, 18, 29], type = "bar", name = "Montreal"),
#             ],
#             layout = (title = "Vamos", barmode="group")
#         )
#     )
# end

run_server(app, "0.0.0.0", debug=true)