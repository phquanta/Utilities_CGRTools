from selenium import webdriver
import time
from icecream import ic
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC

browser = webdriver.Firefox()
browser.get('http://mapper.grzybowskigroup.pl/marvinjs/')
textarea = browser.find_element_by_id("wejscieUzytkownika")
text="[C:1](-[H:2])(-[H:3])(-[O:4]-[H:5])-[C:6](-[H:7])(-[C:8](=[O:9])-[O-:10])-[N+:11](-[H:12])(-[H:13])-[H:14].[O:15](-[H:16])-[H:17].[O:18]=[O:19]>>[O:20](-[H:21])-[O:22]-[H:23].[N+:24](-[H:25])(-[H:26])(-[H:27])-[H:28].[O-:29]-[C:30](=[O:31])-[C:32](=[O:33])-[C:34](-[H:35])(-[H:36])-[O:37]-[H:38]"
textarea.clear()
textarea.send_keys(text)
#textarea.submit()
browser.find_element_by_xpath("//input[@type='submit' and @value='Map reaction!']").click()
#browser.implicitly_wait(24)
time.sleep(4)
textarea.clear()
browser.find_element_by_id("zmapujRysunek").click()
#textarea1 = browser.find_element_by_id("wejscieUzytkownika")
time.sleep(5)
iframe = browser.find_element_by_xpath("//iframe")
browser.switch_to.frame(iframe)
browser.find_element_by_css_selector("[title='Clean 3D']").click()
#WebDriverWait(browser,20).until(EC.element_to_be_clickable((By.XPATH,"//td[@title='Clean 3D')]"))).click()
browser.switch_to.default_content()
textOut=textarea.get_property('value')
ic(textarea)
ic(textOut)
