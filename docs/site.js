document.querySelectorAll("pre > code").forEach((codeBlock) => {
  const pre = codeBlock.parentElement;
  const button = document.createElement("button");
  button.type = "button";
  button.className = "copy-button";
  button.textContent = "Copy";
  button.setAttribute("aria-label", "Copy this R code block");

  button.addEventListener("click", async () => {
    const code = codeBlock.textContent;
    try {
      await navigator.clipboard.writeText(code);
    } catch (_) {
      const textarea = document.createElement("textarea");
      textarea.value = code;
      textarea.style.position = "fixed";
      textarea.style.opacity = "0";
      document.body.appendChild(textarea);
      textarea.select();
      document.execCommand("copy");
      textarea.remove();
    }

    button.textContent = "Copied";
    button.classList.add("is-copied");
    window.setTimeout(() => {
      button.textContent = "Copy";
      button.classList.remove("is-copied");
    }, 1400);
  });

  pre.appendChild(button);
});

const tocLinks = [...document.querySelectorAll(".sidebar a")];
const headings = tocLinks
  .map((link) => document.querySelector(link.getAttribute("href")))
  .filter(Boolean);

if ("IntersectionObserver" in window) {
  const observer = new IntersectionObserver(
    (entries) => {
      const visible = entries.find((entry) => entry.isIntersecting);
      if (!visible) return;
      tocLinks.forEach((link) => {
        link.classList.toggle(
          "is-active",
          link.getAttribute("href") === `#${visible.target.id}`
        );
      });
    },
    { rootMargin: "0px 0px -72% 0px", threshold: 0.05 }
  );
  headings.forEach((heading) => observer.observe(heading));
}
