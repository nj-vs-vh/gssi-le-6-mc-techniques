from pathlib import Path
import re


TEMPLATE_PATH = Path("README.template.md")

template_body = TEMPLATE_PATH.read_text()

INCLUDE_RE = re.compile(
    r"^\!include\[(?P<title>.*?)\]\((?P<filename>.*?)\)(\((?P<programming_lang>.*?)\))?$",
    flags=re.MULTILINE,
)

readme_body = "<!-- DO NOT EDIT, GENERATED AUTOMATICALLY -->\n\n"
current_pos = 0
for include_match in INCLUDE_RE.finditer(template_body):
    readme_body += template_body[current_pos:include_match.start()]
    current_pos = include_match.end()
    print()
    print(f"Found include directive: {include_match.group()}")
    title: str = include_match.group("title")
    filename: str = include_match.group("filename")
    included = Path(filename).read_text()
    programming_lang = include_match.group("programming_lang") or ""
    print(f"Parsed {title = !r} {filename = !r} {programming_lang = !r}")
    included = f"\n\n<details>\n<summary>{title}</summary>\n\n```{programming_lang}\n{included}\n```\n\n</details>\n\n"
    readme_body += included

readme_body += template_body[current_pos:]

Path("README.md").write_text(readme_body)
print("Done")