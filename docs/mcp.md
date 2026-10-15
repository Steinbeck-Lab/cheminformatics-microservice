---
outline: deep
---

# Model Context Protocol (MCP)

## What MCP does

The [Model Context Protocol](https://modelcontextprotocol.io/) gives AI applications
a standard way to discover and call external tools. The Cheminformatics
Microservice exposes its REST operations as MCP tools, allowing an MCP-capable
assistant or agent to use RDKit, CDK, Open Babel, and other chemistry utilities
without implementing a separate HTTP integration for every operation.

MCP is an interface to the microservice, **not an AI model**. The client supplies
the assistant or workflow that chooses which tools to call. The server executes
the chemistry operations and returns their results.

```text
User request
    -> MCP-capable assistant or workflow
    -> MCP tool call
    -> Cheminformatics REST operation
    -> Chemistry toolkit
    -> Result returned to the client
```

For example, a client could turn "Convert ethanol's SMILES to InChI" into a
`generate_inchi` call with `{"smiles": "CCO", "toolkit": "rdkit"}`.
This uses the same implementation as `GET /convert/inchi`.

## Why it is useful

- **Tool discovery:** clients retrieve tool names, descriptions, and JSON input
  schemas instead of maintaining a separate hand-written integration.
- **Natural-language workflows:** assistants can translate requests into
  structured chemistry operations.
- **Multi-step analysis:** a client can retrieve a structure from PubChem,
  standardize it, calculate descriptors, and request a depiction.
- **Shared implementations:** REST and MCP use the same chemistry functions,
  rather than maintaining separate calculation code.
- **Centralized dependencies:** clients do not need their own installation of
  Java, RDKit, CDK, Open Babel, or SURGE.

Example requests for a connected assistant include:

- "Canonicalize this SMILES using RDKit and convert it to InChIKey."
- "Calculate descriptors and identify functional groups for this molecule."
- "Compare these structures using Tanimoto similarity."
- "Detect and remove circular sugars from this structure."
- "Generate structures for the molecular formula C2H6."

The client must support MCP tool use; merely opening the endpoint in a browser
does not provide an assistant or chat interface. Scientific interpretation and
validation of the results remain the responsibility of the user.

## Transport and endpoints

The current lite deployment pins `fastapi-mcp==0.3.4`. It uses **HTTP SSE**, not
Streamable HTTP or stdio.

| Endpoint              | Purpose                                                                        |
| --------------------- | ------------------------------------------------------------------------------ |
| `GET /mcp`            | Persistent SSE connection for MCP clients                                      |
| `POST /mcp/messages/` | Client messages; use the session-specific URL advertised by the SSE connection |
| `GET /docs`           | Interactive REST API documentation                                             |
| `GET /openapi.json`   | REST OpenAPI schema                                                            |
| `GET /health`         | Application health check                                                       |

Configure an SSE-capable MCP client with `https://mcp.naturalproducts.net:8002/mcp`
once DNS and TLS are configured. HTTPS requires a valid certificate and a
TLS-enabled listener on port 8002; DNS mapping alone does not enable it.
Client configuration formats vary; select the SSE transport explicitly where
required. Clients supporting only Streamable HTTP cannot use this adapter as-is.

MCP requests initialize a session, send the initialized notification, and then
use operations such as `tools/list` and `tools/call`. A normal HTTP GET to `/mcp`
opens the event stream; it does not list tools.

## Quick guide: connect a local AI agent

Your agent runs locally, while chemistry calculations run on the MCP server.
Use an agent that supports **remote SSE MCP servers** and a model capable of
tool use. A local model alone is not enough: the agent application must handle
MCP discovery, tool calls, and returning results to the model.

### 1. Check network access

From the computer running your agent:

```bash
curl --fail https://mcp.naturalproducts.net:8002/health
```

The planned hostname is `mcp.naturalproducts.net`. This example will work after
its DNS record and TLS listener are configured and your computer can reach port 8002.

### 2. Add the server to your agent

In an agent with a remote MCP connection form, enter:

| Setting        | Value                            |
| -------------- | -------------------------------- |
| Name           | `cheminformatics`                |
| Transport      | `SSE`                            |
| Server URL     | `https://mcp.naturalproducts.net:8002/mcp` |
| Authentication | None in the current deployment   |

For example, in **VS Code with GitHub Copilot agent mode**, add this workspace
configuration to `.vscode/mcp.json`:

```json
{
  "servers": {
    "cheminformatics": {
      "type": "sse",
      "url": "https://mcp.naturalproducts.net:8002/mcp"
    }
  }
}
```

If the file already exists, add the `cheminformatics` entry to its existing
`servers` object instead of replacing other connections. Start the server
connection using VS Code's MCP controls, review any trust prompt, and enable
its tools in the agent chat tool picker. Model access and subscription
requirements depend on your client.

Other agents use different configuration formats; do not assume this VS Code
JSON can be copied into every client. If an agent supports only stdio or
Streamable HTTP, it needs a compatible bridge or a different server transport.

### 3. Ask the agent to use chemistry tools

Try this prompt:

> Use the cheminformatics MCP tools to canonicalize ethanol's SMILES `CCO`
> and generate its InChI using RDKit. Report the tools used and their results.

The agent should discover and invoke:

```text
generate_canonical_smiles({"smiles": "CCO", "toolkit": "rdkit"})
generate_inchi({"smiles": "CCO", "toolkit": "rdkit"})
```

Expected values are:

```text
Canonical SMILES: CCO
InChI: InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3
```

Approve tool execution if prompted. Check the agent's tool-call trace to confirm
it contacted the server rather than answering from model knowledge.

### 4. Troubleshoot the connection

- **Cannot connect:** check the address, VPN/routing, and inbound TCP port 8002.
  `localhost` means the agent's own computer, not the remote server.
- **Health works but tools do not appear:** confirm SSE support, the `/mcp` URL,
  and `ENABLE_MCP=true` on the server. Review the client's MCP connection logs.
- **Connected but tools are not called:** enable the server's tools, select a
  tool-capable agent/model, and explicitly request tool use.
- **A tool fails:** inspect its discovered input schema and server logs. A
  successful connection does not guarantee external-service or upload support.

::: warning Local agent does not mean local data processing
Molecules, documents, and tool results are sent to the remote chemistry server.
If your agent uses a cloud model, it may also send prompts and results to that
provider. The current HTTP deployment is unencrypted and unauthenticated;
use trusted-network access and avoid sensitive data.
:::

## Exposed tools

The MCP-enabled lite configuration exposes the following **34 tools**. Tool names
come from REST OpenAPI operation IDs; names without explicit operation IDs are
generated by FastAPI. Preserve the exact spelling and capitalization when making
calls. REST paths below are relative to the server base URL and are unversioned
in this configuration.

Tool availability can change with the deployed code and enabled routers. Use
[`tools/list`](#list-tools-and-input-schemas) as the authoritative inventory for
your running server.

### Chemical analysis: 15 tools

| MCP tool                              | REST operation                        | Purpose                                           |
| ------------------------------------- | ------------------------------------- | ------------------------------------------------- |
| `enumerate_stereoisomers`             | `GET /chem/stereoisomers`             | Enumerate possible stereoisomers                  |
| `generate_descriptors`                | `GET /chem/descriptors`               | Calculate molecular descriptors                   |
| `generate_multiple_descriptors`       | `GET /chem/descriptors/multiple`      | Calculate descriptors for multiple molecules      |
| `generate_HOSE_codes`                 | `GET /chem/HOSEcode`                  | Generate HOSE codes                               |
| `standardize_molblock`                | `POST /chem/standardize`              | Standardize a molblock using the ChEMBL pipeline  |
| `check_errors`                        | `GET /chem/errors`                    | Check a SMILES structure for issues               |
| `generate_nplikeness_score`           | `GET /chem/nplikeness/score`          | Calculate natural product likeness                |
| `generate_tanimoto_similarity`        | `GET /chem/tanimoto`                  | Calculate molecular similarity                    |
| `generate_coconut_preprocessing`      | `GET /chem/coconut/pre-processing`    | Generate COCONUT preprocessing information        |
| `generate_classyfire_classifications` | `GET /chem/classyfire/classify`       | Submit a ClassyFire classification                |
| `retrieve_classyfire_results`         | `GET /chem/classyfire/{jobid}/result` | Retrieve classification results by job ID         |
| `filter_molecules`                    | `POST /chem/all_filters`              | Apply selected molecular filters                  |
| `generate_ertl_functional_groups`     | `GET /chem/ertlfunctionalgroup`       | Identify functional groups using Ertl's algorithm |
| `standardize_tautomer`                | `GET /chem/standarizedTautomer`       | Standardize a tautomer                            |
| `retrieve_pubchem_smiles`             | `GET /chem/pubchem/smiles`            | Retrieve SMILES from PubChem by identifier        |

See also the [Chem module](./chem.md).

### Format conversion: 10 tools

| MCP tool                              | REST operation                 | Purpose                                           |
| ------------------------------------- | ------------------------------ | ------------------------------------------------- |
| `generate_2d_coordinates`             | `GET /convert/mol2D`           | Generate 2D molecular coordinates                 |
| `generate_3d_coordinates`             | `GET /convert/mol3D`           | Generate 3D molecular coordinates                 |
| `generate_smiles`                     | `GET /convert/smiles`          | Convert supported input representations to SMILES |
| `generate_canonical_smiles`           | `GET /convert/canonicalsmiles` | Canonicalize SMILES                               |
| `generate_cxsmiles`                   | `GET /convert/cxsmiles`        | Generate CXSMILES                                 |
| `generate_inchi`                      | `GET /convert/inchi`           | Generate InChI                                    |
| `generate_inchikey`                   | `GET /convert/inchikey`        | Generate InChIKey                                 |
| `generate_selfies`                    | `GET /convert/selfies`         | Generate SELFIES                                  |
| `generate_formats`                    | `GET /convert/formats`         | Convert SMILES to multiple molecular formats      |
| `smiles_to_smarts_convert_smarts_get` | `GET /convert/smarts`          | Convert SMILES to SMARTS                          |

See also the [Convert module](./conversions.md).

### Visualization: 2 tools

| MCP tool                           | REST operation   | Purpose                              |
| ---------------------------------- | ---------------- | ------------------------------------ |
| `generate_2d_depiction`            | `GET /depict/2D` | Generate a 2D molecular depiction    |
| `depict_3d_molecule_depict_3D_get` | `GET /depict/3D` | Generate a 3D visualization response |

See also the [Depict module](./depict.md). A client's ability to display a
depiction depends on how it handles the returned content; tool discovery does
not guarantee native image or interactive HTML rendering.

### Structure generation, sugars, and PDFs: 7 tools

| MCP tool                                   | REST operation                      | Purpose                                                  |
| ------------------------------------------ | ----------------------------------- | -------------------------------------------------------- |
| `generate_structures`                      | `GET /tools/generate-structures`    | Generate structures from a molecular formula using SURGE |
| `get_sugar_information`                    | `GET /tools/sugars-info`            | Detect circular and linear sugars                        |
| `remove_linear_sugars`                     | `GET /tools/remove-linear-sugars`   | Remove linear sugars                                     |
| `remove_circular_sugars`                   | `GET /tools/remove-circular-sugars` | Remove circular sugars                                   |
| `remove_linear_and_circular_sugars`        | `GET /tools/remove-sugars`          | Remove both sugar types                                  |
| `extract_pdf_json_tools_extract_json_post` | `POST /tools/extract_json`          | Extract structured content from a PDF                    |
| `extract_pdf_text_tools_extract_text_post` | `POST /tools/extract_text`          | Extract combined text from a PDF                         |

See also the [Tools module](./tools.md).

### Limitations

- **PDF uploads:** PDF operations appear in discovery, but multipart file upload
  through this adapter has not been verified. Use the REST endpoints in `/docs`
  for file uploads until MCP upload support is validated. Docling may need model
  downloads on first conversion.
- **External services:** PubChem and ClassyFire operations require network access
  and depend on the availability and behavior of those services.
- **OCSR:** the lite deployment does not expose DECIMER image-recognition tools.
  A full OCSR-capable image and enabled router are needed.
- **Infrastructure endpoints:** `/health` and `/metrics` are REST endpoints but
  are not part of the 34-tool MCP inventory above.
- **Input details:** use the discovered input schemas for parameter names,
  defaults, and types. Tool registration alone does not prove successful
  execution for every input, especially file uploads and request bodies.

## Execute a tool

Within an initialized SDK session, call a tool by its discovered name:

```python
result = await session.call_tool(
    "generate_canonical_smiles",
    arguments={"smiles": "CCO", "toolkit": "rdkit"},
)
if result.isError:
    raise RuntimeError(f"Tool execution failed: {result}")
for content in result.content:
    if content.type == "text":
        print(content.text)
```

Place this code inside the `ClientSession` block in the previous example, after
`await session.initialize()`. For this input, the text content is the JSON string
`"CCO"`; use `json.loads(content.text)` to decode that particular response.

This verifies a real chemistry operation rather than only a health check.
Different operations can return strings, objects, or other response content;
inspect each tool's description and response before consuming it in a workflow.
