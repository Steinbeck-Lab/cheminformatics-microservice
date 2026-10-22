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

## Quick guide: connect a local AI agent

Your agent runs locally, while chemistry calculations run on the MCP server.
Use an agent that supports **stdio MCP servers** and a model capable of tool use.
The local bridge handles remote SSE connectivity. A local model alone is not
enough: the agent application must handle MCP discovery, tool calls, and returning
results to the model.

### 1. Check network access

From the computer running your agent:

```bash
curl --fail https://mcp.naturalproducts.net:8002/health
```

### 2. Add the server to your agent

For **VS Code with GitHub Copilot agent mode**, configure a local
[`mcp-remote`](https://github.com/geelen/mcp-remote) bridge. VS Code communicates
with it over stdio, and the bridge connects to the server using SSE:

```text
VS Code agent -> local stdio bridge -> remote SSE MCP server
```

The bridge was verified on macOS with initialization, discovery of 34 tools, and
successful RDKit canonical SMILES and InChI calls. The Windows configuration
below follows the same approach but has not been tested on Windows.

Install a supported Node.js LTS release, then verify `node --version` and
`npm --version` in a VS Code terminal. The examples pin `mcp-remote` to the
verified version. `npx -y` downloads and executes that third-party package on
first use, so review and trust it before starting the connection.

Create `.vscode/mcp.json` in your project folder using the configuration for your
operating system below. If the file already exists, add or replace only the
`cheminformatics` entry in its `servers` object, preserving other connections.

**macOS/Linux: `.vscode/mcp.json`**

```json
{
  "servers": {
    "cheminformatics": {
      "type": "stdio",
      "command": "npx",
      "args": [
        "-y",
        "mcp-remote@0.15.1",
        "https://mcp.naturalproducts.net:8002/mcp",
        "--transport",
        "sse-only"
      ]
    }
  }
}
```

If VS Code cannot find `npx`, run `command -v npx` in a terminal and use the
resulting absolute path for `command`. For example, Apple Silicon Homebrew
commonly installs it at `/opt/homebrew/bin/npx`. Restart VS Code after installing
Node.js so it inherits the updated environment.

**Windows: `.vscode/mcp.json`**

Use `cmd.exe` to launch the Windows `npx` command:

```json
{
  "servers": {
    "cheminformatics": {
      "type": "stdio",
      "command": "cmd.exe",
      "args": [
        "/d",
        "/c",
        "npx",
        "-y",
        "mcp-remote@0.15.1",
        "https://mcp.naturalproducts.net:8002/mcp",
        "--transport",
        "sse-only"
      ]
    }
  }
}
```

#### Start the connection and enable Copilot tools

1. Open the project folder containing `.vscode/mcp.json` and save the file.
2. Run **MCP: List Servers** from the Command Palette, select
   **cheminformatics**, and choose **Start Server** or **Restart Server**.
   The editor's inline Start button is not required.
3. Review any trust prompt. If the old connection remains cached, run
   **Developer: Reload Window** and start the connection again.
4. Open GitHub Copilot Chat, select **Agent** mode, and enable the
   **cheminformatics** tools through **Configure Tools**.

Model access and subscription requirements depend on your client. Other agents use different configuration formats; do not assume this VS Code JSON can be copied into every client. Configure the same stdio command and arguments in the format your agent supports.

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
